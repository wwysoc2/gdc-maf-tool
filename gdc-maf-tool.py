import os
import re
import sys
import json
import csv
import requests
import hashlib
import gzip
import argparse
import datetime
import logging

def arg_parse():
    parser = argparse.ArgumentParser(
        description = '----GDC MAF Concatenation Tool v1.0----',
              usage = 'python gdc-maf-cat.py <-m MANIFEST or -p PROJECT_ID>')
    parser.add_argument('-m', '--manifest', action = "store",
        help = 'Specify MAF files with GDC Manifest')
    parser.add_argument('-p', '--project', action = "store",
        help = 'Specify MAF files by project')
    parser.add_argument('-o', '--output', action = "store", 
        dest = 'o', type = str, default = "outfile.maf",
        help = 'Designates a name for the output file')
    parser.add_argument('--metrics_only', action = "store_true",
        dest = 'mo', help = 'Only return metrics, do not output concatenated MAF.')
    args = parser.parse_args()
    return args

def main(args):
    '''
    Retrieves and parses the arguments
    '''
    global use_manifest, output_file, manifest_path, project_string, mo
    if args.manifest: 
        use_manifest = True
        manifest_path = args.manifest
    if args.project: 
        use_manifest = False
        project_string = args.project
    if args.o: output_file = args.o
    if args.manifest and args.project:
        error_parse("both_argue")
    if not args.manifest and not args.project:
        error_parse("no_argue")
    mo = False
    if args.mo:
        mo = True
        if args.mo and args.o:
            info_parse("Warning: Metrics-Only mode selected. No data will be written to {}".format(args.o))

def error_parse(code):
    '''
    Generates the ERROR messages
    '''
    error = {
        "bad_manifest": "Input must be valid GDC Manifest. " \
        "\tGo to https://portal.gdc.cancer.gov/ to download a manifest",
        "no_result": "Query produced no results",
        "no_argue": "No argument detected, please use the -p or -m flags",
        "both_argue": "Must choose either -p OR -m, not both.",
        "md5sum_mis": "Expected md5sum does not match file's md5sum value",
        "max_retry" : "Maximum retries exceeded" 
    }
    logger = get_logger("GDC-MAF-Tool")
    logger.error(error[code])
    sys.exit(2)

def info_parse(message):
    '''
    Generates the INFO messages
    '''
    logger = get_logger("GDC-MAF-Tool")
    logger.info(message)    


def strip_maf_header(maf_file):
    '''
    Removes the MAF header
    '''
    maf_list = []
    header = {}
    for line in maf_file:
        if line[0] != "#":
            maf_list.append(line)    
        elif line[0] == "#":
            key, value = line[1:].strip().split(" ")
            header[key] = value

    return maf_list,header

def merge_maf_headers(header_json, output_file):
    '''
    Merges the headers from all MAFs to cut down on duplication
    '''
    keylist = []
    for single in header_json:
        keylist.extend(single.keys())
        keylist = list(set(keylist))
    o = open(output_file, "w")
    for key in keylist:
        one_key = []
        for single in header_json:
            one_key.append(single[key])
            one_key = list(set(one_key))
        header_line = "#{} {}".format(key, ";".join(one_key))
        o.write(header_line + "\n")
    o.close

def jsonify_maf(maf_file):
    '''
    Converts MAF TSV to dict, requires header is stripped.
    '''
    master_dict = []
    keys = maf_file[0].strip().split("\t")
    for line in maf_file[1:]:
        split_line = line.strip().split("\t")
        one_line_dict = dict(zip(keys, split_line))
        master_dict.append(one_line_dict)
    return master_dict, keys

def back_to_tsv(full_dict, col_order, prefix):
    '''
    Converts full concatenated dict to TSV for writing out
    '''
    dict_writer = csv.DictWriter(open("{}".format(prefix), "a"), col_order, delimiter='\t')            
    dict_writer.writeheader()
    dict_writer.writerows(full_dict)
    info_parse("Concatenated MAF written to: {}".format(prefix))
    
def read_in_manifest(manifest_path):
    '''
    Reads in a GDC Manifest to parse out UUIDs
    '''
    manifest_file = open(manifest_path, "r").read().splitlines()
    id_list = []
    if manifest_file[0].strip().split("\t")[0] != "id":
        error_parse("bad_manifest")
    for line in manifest_file[1:]:
        id_list.append(line.strip().split("\t")[0])
    return id_list

def retrieve_ids_by_project(provided, project):
    '''
    Retrieves IDs when provided a project_id or list of UUIDs
    '''
    id_list = []
    endpt = "https://api.gdc.cancer.gov/files"
    filters = [
                ("files.data_format",["MAF"]),
                ("files.data_type",["Masked Somatic Mutation"])]
    if project == True:
        filters.append(("cases.project.project_id", provided.split(",")))
    else:
        filters.append(("files.file_id", provided))
    filters_gdc = {"op":"and", "content":[]}
    for field, value in filters:
        filt_core = {"field": field, "value": value}
        single_filt = {"op": "in", "content": filt_core}
        filters_gdc["content"].append(single_filt)
    params = {
        "filters": json.dumps(filters_gdc),
        "fields" : "file_id,md5sum,file_name",
        "format" : "JSON",
        "size"   : "10000"
    }  
    response = requests.get(endpt, params= params)
    out_hits = json.loads(response.content)["data"]["hits"]
    if len(out_hits) == 0:
        error_parse("no_result")
    for file_entry in out_hits:
        single_dict = dict(zip(["file_id", "md5sum", "file_name"],
            [file_entry["file_id"], file_entry["md5sum"], file_entry["file_name"]]))
        id_list.append(single_dict)
    return id_list

def download_maf(single_maf_dict, tmpdir):
    '''
    Downloads each MAF file and stores in tmp directory
    '''
    file_id, exp_md5 = single_maf_dict["file_id"], single_maf_dict["md5sum"]
    retry = True
    retry_num = 0
    while retry == True and retry_num < 3:
        data_endpt = "https://api.gdc.cancer.gov/data/{}".format(file_id)
        info_parse("> {} | Downloading File | {} |".format(datetime.datetime.now(), file_id))
        response = requests.get(data_endpt, headers = {"Content-Type": "application/json"})
        if response.status_code == 200:
            retry = False
        else:
            retry_num += 1
            info_parse("> -- Retrying Download...")
    if retry == False:
        response_head_cd = response.headers["Content-Disposition"]
        file_name = re.findall("filename=(.+)", response_head_cd)[0]
        with open("/".join([tmpdir, file_name]), "wb") as output_file:
            output_file.write(response.content)
        check_md5sum(file_name, exp_md5, tmpdir)
    elif retry_num == 3:
        error_parse("max_retry")

def download_run(id_list):
    '''
    Runs MAF download for multiple and performed per-session tasks
    '''
    tmpdir = "tmpMAF_" + str(datetime.datetime.now()).split(" ")[0]
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    for single_maf in id_list:
        download_maf(single_maf, tmpdir)
    info_parse("> All MAF Downloads Complete")
    return id_list, tmpdir

def check_md5sum(file_name, exp_md5, tmpdir):
    '''
    Checks the MD5SUM matches the one in the GDC index
    '''
    hash_md5 = hashlib.md5()
    with open("/".join([tmpdir, file_name]), "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    if exp_md5 != hash_md5.hexdigest():
        error_parse("md5sum_mis")

def output_json_file(json_obj, name):
    json.dump(json_obj, open("{}.json".format(name), "w"), indent=4)

def calc_basic_metrics(maf_dict):
    gene_list = []
    case_list = []
    mutation_list = []
    for entry in maf_dict:
        gene_list.append(str(entry["Entrez_Gene_Id"]))
        case_list.append(entry["case_id"])
        mutation = ":".join([str(entry["Entrez_Gene_Id"]),str(entry["Start_Position"]),str(entry["End_Position"]),str(entry["Chromosome"])])
        mutation_list.append(mutation)
    for message in ["> ---Basic-Metrics---",
        "> {} Genes".format(str(len(set(gene_list)))),
        "> {} Cases".format(str(len(set(case_list)))),
        "> {} Mutations".format(str(len(set(mutation_list))))]:
        info_parse(message)

def date_file_format(datetime_str):
    new_str = str(datetime_str).split('.')[0].replace(":", "_").replace(" ", "_")
    return new_str

def get_logger(name):
    global loggers
    if loggers.get(name):
        return loggers.get(name)
    else:
        logger = logging.getLogger(name)
        handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter('[%(asctime)s][%(name)10s][%(levelname)7s] %(message)s')
        handler.setFormatter(formatter)
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)
        loggers[name] = logger
        return logger

def execute():
    main(arg_parse())
    cat_maf = []
    if use_manifest == True:
        maf_ids_only = read_in_manifest(manifest_path) 
        maf_ids = retrieve_ids_by_project(maf_ids_only, False)
    else:
        maf_ids =  retrieve_ids_by_project(project_string, True)

    id_list, tmpdir  = download_run(maf_ids)

    header = []
    for single_maf in id_list:
        maf_list, header_json = strip_maf_header(gzip.open("/".join([tmpdir,single_maf["file_name"]]), "r"))
        header.append(header_json)
        jsonified, keys = jsonify_maf(maf_list)
        cat_maf += jsonified

    calc_basic_metrics(cat_maf)
    
    if mo == False:
        merge_maf_headers(header, output_file)
        back_to_tsv(cat_maf, keys, output_file)

loggers = {}
execute()

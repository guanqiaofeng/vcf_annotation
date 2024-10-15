#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 Copyright (c) 2024, Ontario Institute for Cancer Research (OICR).

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as published
 by the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program. If not, see <https://www.gnu.org/licenses/>.

 Author: Guanqiao Feng <gfeng@oicr.on.ca>
 """

from datetime import date
import argparse
import copy
import hashlib
import json
import os
import sys
import uuid
import yaml


def calculate_size(file_path):
    return os.stat(file_path).st_size


def calculate_md5(file_path):
    md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            md5.update(chunk)
    return md5.hexdigest()


def rename_file(f, payload, seq_experiment_analysis_dict, date_str):
    experimental_strategy = payload['experiment']['experimental_strategy'].lower()

    if f.endswith('.vcf.gz'):
        file_ext = 'vcf.gz'
    elif f.endswith('.vcf.gz.tbi'):
        file_ext = 'vcf.gz.tbi'
    else:
        sys.exit('Error: unknown aligned seq extention: %s' % f)

    variant_type = ''
    if 'snv' in seq_experiment_analysis_dict['files'][0]['fileName'].lower():
        variant_type = 'snv'
    elif 'indel' in seq_experiment_analysis_dict['files'][0]['fileName'].lower():
        variant_type = 'indel'
    elif 'rearrangement' in seq_experiment_analysis_dict['files'][0]['fileName'].lower():
        variant_type = 'sv'
    elif 'copy_number' in seq_experiment_analysis_dict['files'][0]['fileName'].lower():
        variant_type = 'cnv'
    else:
        sys.exit('Error: unknown variant type: %s' % f)

    new_name = "%s.%s.%s.%s.%s.%s.%s.%s.%s.%s" % (
        payload['studyId'],
        seq_experiment_analysis_dict['samples'][0]['specimen']['donorId'],
        seq_experiment_analysis_dict['samples'][0]['sampleId'],
        experimental_strategy,
        date_str,
        seq_experiment_analysis_dict['workflow']['workflow_short_name'],
        seq_experiment_analysis_dict['variant_class'].replace(",", "-").lower(),
        variant_type,
        'annotated',
        file_ext
    )

    new_dir = 'out'
    try:
        os.mkdir(new_dir)
    except FileExistsError:
        pass

    dst = os.path.join(os.getcwd(), new_dir, new_name)
    os.symlink(os.path.abspath(f), dst)

    return dst


def get_files_info(file_to_upload):

    if file_to_upload.endswith('.tbi'):
        data_type = "VCF Index"
        # Check for variations in the third-last part (e.g., .snv., .indel., etc.)
        if '.snv.' in file_to_upload:
            data_category = 'Simple Nucleotide Variation'
        elif '.indel.' in file_to_upload:
            data_category = 'Simple Nucleotide Variation'
        elif '.sv.' in file_to_upload:
            data_category = 'Structural Variation'
        elif '.cnv.' in file_to_upload:
            data_category = 'Copy Number Variation'
        else:
            raise ValueError(f"Data type not recognized for file: {file_to_upload}")
    else:
        # Non-index file (.vcf.gz)
        if '.snv.' in file_to_upload:
            data_type = 'Raw SNV Calls'
            data_category = 'Simple Nucleotide Variation'
        elif '.indel.' in file_to_upload:
            data_type = 'Raw InDel Calls'
            data_category = 'Simple Nucleotide Variation'
        elif '.sv.' in file_to_upload:
            data_type = 'Raw SV Calls'
            data_category = 'Structural Variation'
        elif '.cnv.' in file_to_upload:
            data_type = 'Raw CNV Calls'
            data_category = 'Copy Number Variation'
        else:
            raise ValueError(f"Data type not recognized for file: {file_to_upload}")

    return {
        'fileName': os.path.basename(file_to_upload),
        'fileType': 'VCF' if file_to_upload.split(".")[-2] == 'vcf' else 'TBI',
        'fileSize': calculate_size(file_to_upload),
        'fileMd5sum': calculate_md5(file_to_upload),
        'fileAccess': 'controlled',
        'dataType': data_type, # update
        'info': {
            'data_category': data_category # update
            }
    }


def get_sample_info(sample_list):
    samples = copy.deepcopy(sample_list)
    for sample in samples:
        for item in ['info', 'sampleId', 'specimenId', 'donorId', 'studyId']:
            sample.pop(item, None)
            sample['specimen'].pop(item, None)
            sample['donor'].pop(item, None)

    return samples


def main(args):
    with open(args.seq_experiment_analysis, 'r') as f:
        seq_experiment_analysis_dict = json.load(f)

    pipeline_info = {}
    if args.pipeline_yml:
        with open(args.pipeline_yml, 'r') as f:
            pipeline_info = yaml.safe_load(f)

    updated_pipeline_info = {}
    for key, value in pipeline_info.items():
        new_key = key.split(":")[-1]
        updated_pipeline_info[new_key] = value

    for key, value in updated_pipeline_info.items():
        for sub_key, sub_value in value.items():
            value[sub_key] = str(sub_value)
        updated_pipeline_info[key] = value

    payload = {
        'analysisType': {
            'name': seq_experiment_analysis_dict['analysisType']['name']
        },
        'studyId': seq_experiment_analysis_dict.get('studyId'),
        'info': {},
        'workflow': {
            'workflow_name': args.wf_name,
            'workflow_version': args.wf_version,
            'genome_build': args.genome_build,
            'run_id': args.wf_run,
            'session_id': args.wf_session,
            'pipeline_info': updated_pipeline_info,
            'inputs': [
                {
                    'analysis_type': 'sequencing_experiment',
                    'input_analysis_id': seq_experiment_analysis_dict.get('analysisId')
                }
            ]
        },
        'variant_class': seq_experiment_analysis_dict.get('variant_class'),
        'files': [],
        'samples': get_sample_info(seq_experiment_analysis_dict.get('samples')),
        'experiment': {},
    }
    if args.genome_annotation:
        payload['workflow']['genome_annotation'] = args.genome_annotation

    # pass `info` dict from seq_experiment payload to new payload
    if 'info' in seq_experiment_analysis_dict and isinstance(seq_experiment_analysis_dict['info'], dict):
        payload['info'] = seq_experiment_analysis_dict['info']
    else:
        payload.pop('info')

    payload['experiment'].update(seq_experiment_analysis_dict.get('experiment', {}))

    if 'library_strategy' in payload['experiment']:
        experimental_strategy = payload['experiment'].pop('library_strategy')
        payload['experiment']['experimental_strategy'] = experimental_strategy

    # get file of the payload
    date_str = date.today().strftime("%Y%m%d")
    for f in args.files_to_upload:
        renamed_file = rename_file(f, payload, seq_experiment_analysis_dict, date_str)
        payload['files'].append(get_files_info(renamed_file))

    with open("%s.%s.payload.json" % (str(uuid.uuid4()), args.wf_name.replace(" ","_")), 'w') as f:
        f.write(json.dumps(payload, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tool: payload-rnaseq-alignment')
    parser.add_argument("-f", "--files_to_upload", dest="files_to_upload", type=str, required=True,
                        nargs="+", help="Aligned reads files to upload")
    parser.add_argument("-a", "--seq_experiment_analysis", dest="seq_experiment_analysis", required=True,
                        help="Input analysis for sequencing experiment", type=str)
    parser.add_argument("-w", "--wf_name", dest="wf_name", required=True, help="Workflow name")
    parser.add_argument("-v", "--wf_version", dest="wf_version", required=True, help="Workflow version")
    parser.add_argument("-r", "--wf_run", dest="wf_run", required=True, help="workflow run ID")
    parser.add_argument("-s", "--wf_session", dest="wf_session", required=True, help="workflow session ID")
    parser.add_argument("-b", "--genome_build", dest="genome_build", help="Genome build")
    parser.add_argument("-n", "--genome_annotation", dest="genome_annotation", help="Genome annotation")
    parser.add_argument("-p", "--pipeline_yml", dest="pipeline_yml", required=False, help="Pipeline info in yaml")

    args = parser.parse_args()

    main(args)

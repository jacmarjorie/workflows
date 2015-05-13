#!/usr/bin/env python

import sys
import subprocess
from bioblend import galaxy


def buildParams(read_group_dict_tumor, read_group_dict_normal):

    param_set = {}
    
    param_set['327'] = {'params|readGroup|rgid':read_group_dict_tumor['rgid'], 'params|readGroup|rgpu':read_group_dict_tumor['rgpu'], \
                            'params|readGroup|rglb':read_group_dict_tumor['rglb']}
    param_set['328'] = {'params|readGroup|rgid':read_group_dict_normal['rgid'], 'params|readGroup|rgpu':read_group_dict_normal['rgpu'], \
                            'params|readGroup|rglb':read_group_dict_normal['rglb']}

    return param_set


def buildDatamap():

    datamap = {}

### For small workflow
# datamap['345'] = {'src':'ldda', 'id':'28053eaab7322e74'}
# datamap['346'] = {'src':'ldda', 'id':'0495671efad17eb0'}
# datamap['347'] = {'src':'ldda', 'id':'96586493003fa534'}

### Reference genome
    datamap['351'] = {'src':'ldda', 'id':'28053eaab7322e74'}
### Tumor R1
    datamap['349'] = {'src':'ldda', 'id':'0495671efad17eb0'}
### Tumor R2
    datamap['350'] = {'src':'ldda', 'id':'96586493003fa534'}
### Normal R1
    datamap['352'] = {'src':'ldda', 'id':'35c5e7060e4be4ca'}
### Normal R2
    datamap['353'] = {'src':'ldda', 'id':'0e2f009200238f96'}
### COSMIC
    datamap['358'] = {'src':'ldda', 'id':'989b1cf8cc875a83'}
### Centromere
    datamap['359'] = {'src':'ldda', 'id':'3e4cf71537771c57'}
### DBSnp
    datamap['356'] = {'src':'ldda', 'id':'990f054d8be549c8'}
### Interval List
    datamap['357'] = {'src':'ldda', 'id':'d445d9dacb99fa13'}
### Mills
    datamap['354'] = {'src':'ldda', 'id':'dbe26252350d9c9f'}
### 1000G
    datamap['355'] = {'src':'ldda', 'id':'3f312d6d70e5d08f'}

    return datamap

def getFastqId(fastq):

    read_groups = []
    
    if fastq.lower().endswith('.gz'):
        cmd = "gzip -cd " + fastq + " | head -1"
    elif fastq.lower().endswith('.fastq') or fastq.lower().endswith('.fq'):
        cmd = "head -1 " + fastq
    else:
        raise Exception("File must end in fastq.gz, fastq, or fq.")

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    while True:
        line = p.stdout.readline()
        if line != '':
            read_groups.extend(line.rstrip('\n').split(':'))
        else:
            break

    return read_groups

def parseFastqId(fastq_id):
    
    read_group_dict = {}

    read_group_dict['rgid'] = ':'.join(fastq_id[:4])
    read_group_dict['rglb'] = ':'.join(fastq_id[:4])
    read_group_dict['rgpu'] = fastq_id[9]

    return read_group_dict


def main():

### API Key from Galaxy.
    api_key = "29ac006cf3af208b5ac1e18800d7894a"

### Find workflow_id.
#print(gi.workflows.get_workflows())
    workflow_id = "dbe26252350d9c9f"
### Workflow_id for smaller test workflow.
#workflow_id = "b349f61fb5a72a65"

### Galaxy Instance URL
    galaxy_url = "http://172.20.15.50:3837/"

    fastq_id_tumor = getFastqId(sys.argv[1])
    fastq_id_normal = getFastqId(sys.argv[2])
    rg_dict_tumor = parseFastqId(fastq_id_tumor)
    rg_dict_normal = parseFastqId(fastq_id_normal)
    datamap = buildDatamap()
    param_set = buildParams(rg_dict_tumor, rg_dict_normal)

    gi = galaxy.GalaxyInstance(url=galaxy_url, key=api_key)

### Find tool id's based on workflow_id.
    print("Inputs to the workflow include: ")
    for key in gi.workflows.show_workflow(workflow_id)['inputs']:
        print(key + '\t' + gi.workflows.show_workflow(workflow_id)['inputs'][key]['label'])


    gi.workflows.run_workflow(workflow_id, datamap, params=param_set, history_name=rg_dict_tumor['rgid'], import_inputs_to_history=False)

if __name__ == "__main__":
    main()

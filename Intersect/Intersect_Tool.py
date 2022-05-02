#!/usr/bin/env python
# coding: utf-8

# # Intersect Tool

# In[1]:


from subprocess import run
import os
import sys
import argparse

# In[2]:


def intersectTool(bam_path, bed_path, intersect=False):
    """
    param str bam_path: <filename>.bam
    param str bed_path: <filename>.bed
    param bool intersect: True if you want reads only intersecting the bed
                          False if you want reads not intersecting the bed
                          default = False, so entering in two arguments will assume 
                                    you want no intersects
    
    """
    if intersect == True:
        output_bam = f"{os.path.splitext(bam_path)[0]}.intersected.bam"
        olap_tool = f"bedtools intersect -a {bam_path} -b {bed_path} > {output_bam}"
        run(olap_tool, shell=True)
    else:
        output_bam = f"{os.path.splitext(bam_path)[0]}.nonintersected.bam"
        olap_tool = f"bedtools intersect -a {bam_path} -b {bed_path} -v > {output_bam}"
        run(olap_tool, shell=True)
    return(output_bam)


# In[3]:


def indexFile(bam_output):
    """
    param str bam_output: copy/paste output_bam from intersectTool()
    
    """
    sammy_tool = f"samtools index {bam_output} {bam_output}.bai"
    run(sammy_tool, shell=True)
    
def parse_cmdline_params():
    """
    parses command line arguments
    
    """
    parse = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parse.add_argument("-ba", "--bam_input",
                       help = "input bam file",
                       type = str,
                       required = True)
    parse.add_argument("-be", "--bed_input",
                       help = "input bed file",
                       type = str,
                       required = True)
    parse.add_argument("-bo", "--boolean_input",
                       help = "true for 3 args, false for 2",
                       required = False,
                       action = "store_true")
    
    opts, unknown = parse.parse_known_args(args=sys.argv[1:])

    return opts
    
if __name__ == "__main__":
    opts = parse_cmdline_params()
    sammy = intersectTool(opts.bam_input, opts.bed_input, opts.boolean_input)
    indexFile(sammy)

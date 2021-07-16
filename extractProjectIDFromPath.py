#!/usr/bin/env python3

import sys

cwd=sys.argv[1].split("/")
projectNo=[x for x in cwd if x.find("Proj_")>-1]
if(len(projectNo)>0):
    print(projectNo[-1])
else:
    print(cwd[-1])

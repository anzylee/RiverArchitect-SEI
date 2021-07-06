import os, sys

import shutil

site = 'sfe_209'
sites = [322,95,82,4523,81,24,25,316,2248,221,209]
ver = [1,1,3,1,0,1,0,9,1,1,2]
ver_ind = sites.index(int(site.split('_')[1]))
case_of_interest = [site, site+'_RB_v'+str(ver[ver_ind])]

original_dir = 'F:/RiverArchitect-development/SHArC/SHArea/' + case_of_interest[1]
target_dir = 'F:/RiverArchitect-development/SHArC/SHArea/' + case_of_interest[0]

for original in os.listdir(original_dir):
    target = original.replace(case_of_interest[1], case_of_interest[0]+'_RB')
    shutil.copyfile(os.path.join(original_dir,original),
                    os.path.join(target_dir,target))

print(site + ' - done')
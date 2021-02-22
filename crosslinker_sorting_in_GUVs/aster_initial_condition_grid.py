#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 15:22:10 2020

@author: redford
"""

"""
If you're reading this, I'm sorry this isn't better organized. This is the script
to create initial configurations and write a file to submit to slurm to run AFINES
simuations for Bashirzadeh et al. 
"""



import numpy as np
from matplotlib import pyplot as plt
import random

storage_root = '/project2/gardel/steven/AFINES/inputs/'
hh = 0

for alpha_conc in range(3,10,6):
    for circle_radius in range(6,13,2):
        
        hh = hh + 1
        #print(hh)
        #print(on_rte)
        #print('alpha = ', alpha_conc)
        #print('nbeads = ', nbeads)
        
        #alpha_conc = 0
        fasc_conc = 3
        r = 11
        area = np.pi*np.square(r)
        nfascin = int(np.round(fasc_conc*area))
        nalpha = int(np.round(alpha_conc*area))
        rfil = 15
        r_rand = 3
        link_len = 1
        nbeads = 31
        fascin_len = 0.06
        alpha_len = 0.36
        nrep = 2
        nfil = 50
        #off_rte = 0.01
        #circle_radius = 12
        #alignment = 0.1*algn
        #a_m_stf = 1
        #a_m_bind = 5
        #am_on = 15*on
        #a_m_bend = 0.006*bnd
        
        
        
        init_x = np.linspace(-rfil,rfil,nbeads)
        init_y = np.zeros(nbeads)
        
        #init_x = 10
        #init_y = 0
        
        #plt.plot(init_x,init_y, 'k-o')
        #plt.show()
        
        xes = np.zeros((nfil,nbeads))
        yes = np.zeros((nfil,nbeads))
        inds = np.zeros((nfil,nbeads))
        
        fascin = np.zeros((2,nfascin))
        alpha = np.zeros((2,nalpha))
        
        
        act_tail = 'actin_%d.txt' %hh
        path = storage_root + act_tail
        f = open(path, 'w+')   
        for ii in range(nfil):
            ang = np.around(random.uniform(0,2*np.pi),4)
            #fil_ang = np.mod((ang + np.pi),2*np.pi)
            #print(fil_ang)
            xes[ii,:] = init_x*np.cos(ang) + init_y*np.sin(ang) + random.uniform(-r_rand,r_rand)
            yes[ii,:] = -init_x*np.sin(ang) + init_y*np.cos(ang) + random.uniform(-r_rand,r_rand)
            inds[ii,:] = int(ii)    
            
            for tt in range(nbeads):
                
                x = np.around(xes[ii,tt],3)
                y = np.around(yes[ii,tt],3)
                idd = inds[ii,tt]
            
                line = '%s\t' %x +  '%s\t' %y + '0.5\t' + '%s\n' %idd
            
                f.write(line)
        f.close()
        
                
                
                
        """      
        for jj in range(nfil):
            plt.plot(xes[jj,:],yes[jj,:],'k-o')
        plt.xlim([-r,r])
        plt.ylim([-r,r])
        plt.axis('equal')
        plt.show()
        """
        
        fasc_tail = 'fascin_%d.txt' %nfascin
        path = storage_root + fasc_tail
        f = open(path, 'w+')          
        for ll in range(nfascin):
            x = np.around(random.uniform(-r,r),3)
            ymax = np.sqrt(np.square(r) - np.square(x))
            y = np.around(random.uniform(-ymax,ymax),3)
            ang = random.uniform(0,2*np.pi)
            dx = np.around(np.cos(ang)*fascin_len,3)
            dy = np.around(np.sin(ang)*fascin_len,3)
            #if ll%100 == 0:
                #print(dx)
                #print(dy)
            
            
            fascin[0,ll] = x
            fascin[1,ll] = y
            
            line = '%s\t' %x +  '%s\t' %y + '%s\t' %dx + '%s\t' %dy + '-1\t' + '-1\t' + '-1\t' + '-1\n'
            
            f.write(line)
        f.close()
        
        
        
        alph_tail = 'alpha_%d.txt' %nalpha
        path = storage_root + alph_tail
        f = open(path, 'w+')  
        for mm in range(nalpha):
            x = random.uniform(-r,r)
            ymax = np.sqrt(np.square(r) - np.square(x))
            y = random.uniform(-ymax,ymax)
            ang = random.uniform(0,2*np.pi)
            dx = np.around(np.cos(ang)*alpha_len,3)
            dy = np.around(np.sin(ang)*alpha_len,3)
            #if mm%100 == 0:
                #print(dx)
                #print(dy)
            alpha[0,mm] = x
            alpha[1,mm] = y
        
            line = '%s\t' %x +  '%s\t' %y + '%s\t' %dx + '%s\t' %dy + '-1\t' + '-1\t' + '-1\t' + '-1\n'
            
            f.write(line)
        f.close()
        
        """
        plt.scatter(fascin[0,:],fascin[1,:], color = 'r')
        plt.scatter(alpha[0,:], alpha[1,:], color = 'c')
        plt.xlim([-r,r])
        plt.ylim([-r,r])
        plt.axis('equal')
        plt.show()
        """
        
        
        
        
        
        path = '/project2/gardel/steven/AFINES/0930_aster_sort_fasc_only_bnd_and_on_%d.csh' % hh
        date = '093020'
        dir_1 = '"/project2/gardel/steven/data/out_%s' % date
        dir_2 = dir_1 + '/aster_sort_fasc_only_bnd_and_on'
        dir_3 = dir_2 + '_cond_%s' % hh
        the_dir = dir_3 + '_trial_$SLURM_ARRAY_TASK_ID"'
            
        actin_path = 'inputs/%s' %act_tail
        fasc_path = 'inputs/%s' %fasc_tail
        alph_path = 'inputs/%s' %alph_tail
        
        
        l1 = '#!/bin/sh \n'
        l2 = '#SBATCH --job-name=sort%d \n' % nbeads
        l3 = '#SBATCH --output=' + date + '_aster_sort_full_long_%d.out \n' % hh
        l4 = '#SBATCH --account=pi-gardel \n'
        l5 = '#SBATCH --partition=broadwl \n'
        l6 = '#SBATCH --qos=normal \n'
        l7 = '#SBATCH --ntasks=1 \n'
        l8 = '#SBATCH --array=1-%s \n' % nrep
        l9 = '#SBATCH --time=36:00:00 \n'
        l10 = '\n'
        l11 = 'direc=%s \n' % the_dir
        l12 = 'rst_tail="/data/config_full.cfg" \n'
        l13 = 'mkdir -p $direc \n'
        l14 = '\n'
        l15 = '/project2/dinner/afines-integrate -c config/sorting_aster_0924.cfg --dir $direc --myseed $RANDOM --actin_in %s' %actin_path
        l15 = l15 + ' --a_motor_in %s' %fasc_path
        #l15 = l15 + ' --a_m_kalign %s' %alignment
        #l15 = l15 + ' --p_m_koff %s' %off_rte
        #l15 = l15 + ' --a_m_koff %s' %am_on
        l15 = l15 + ' --a_motor_density %s' %fasc_conc
        l15 = l15 + ' --p_motor_density %s' %alpha_conc
        #l15 = l15 + ' --a_motor_stiffness %s' %a_m_stf
        #l15 = l15 + ' --p_m_bend %s' %bend
        #l15 = l15 + ' --a_m_bend %s' %a_m_bend
        l15 = l15 + ' --nmonomer %s' %nbeads
        l15 = l15 + ' --npolymer %s' %nfil
        l15 = l15 + ' --circle_radius %s' %circle_radius
        l15 = l15 + ' --p_motor_in %s \n' %alph_path
        
        """
        f = open(path, 'w+')
        
        for ll in range(15):
            line = 'l%d' % (ll + 1)
            f.write(eval(line))
        
        f.close()
        """

print("Done! :D")


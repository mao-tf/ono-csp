#!/bin/bash 

g16 < test_t.gjf > test_t.log 

g16 < test_t_m1.gjf > test_t_m1.log 

g16 < test_t_m2.gjf > test_t_m2.log 

g16 < test_p.gjf > test_p.log 

g16 < test_p_m1.gjf > test_p_m1.log 

g16 < test_p_m2.gjf > test_p_m2.log 

formchk test.chk test.fch 


#sleep 5 

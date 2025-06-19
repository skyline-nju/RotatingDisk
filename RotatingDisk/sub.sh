#! /bin/bash

export PATH=/home/sli/yduan/local/gcc-14.2.0/bin:$PATH
export LD_LIBRARY_PATH=/home/sli/yduan/local/gcc-14.2.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/sli/yduan/local/gcc-14.2.0/lib64:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/home/sli/yduan/local/gcc-14.2.0/include:$C_INCLUDE_PATH

T=0.3
sigma=0
J1=1
h=0.1
n_step=2900000
snap_dt=10000
IC=resume
seed=3000
for omega0 in 0
do
	job_name="NRKM1024_T${T}_o${omega0}_S${seed}"
	cat <<EOF > ${job_name}.pbs
#!/bin/bash
#PBS -N ${job_name}
#PBS -l nodes=1:ppn=1
#PBS -j oe


cd \$PBS_O_WORKDIR

./a.out $T $sigma $omega0 $J1 $h $seed $snap_dt $n_step $IC
EOF
sleep 0.25
qsub ${job_name}.pbs
done

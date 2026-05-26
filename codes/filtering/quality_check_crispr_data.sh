#!/bin/bash

# ###### Zona de Parámetros de solicitud de recursos a SLURM ############################
#
#SBATCH --job-name=fastqc_crispr	#Nombre del job
#SBATCH -p short			#Cola a usar, Default=short (Ver colas y límites en /hpcfs/shared/README/partitions.txt)
#SBATCH -N 1				#Nodos requeridos, Default=1
#SBATCH -n 1				#Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=4		#Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --mem=12000		#Memoria en Mb por CPU, Default=2048
#SBATCH --time=40:00:00			#Tiempo máximo de corrida, Default=2 horas
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL			
#SBATCH -o fastqc_crispr_trimmed.o%j		#Nombre de archivo de salida
#
########################################################################################

# ################## Zona Carga de Módulos ############################################


########################################################################################


# ###### Zona de Ejecución de código y comandos a ejecutar secuencialmente #############
sleep 60
host=`/bin/hostname`
date=`/bin/date`
echo "Soy un JOB de descarga"
echo "Corri en la maquina: "$host
echo "Corri el: "$date
echo -e "Ejecutar Script de python para analizar los CTLs de ZIBO: +18K predios\n"

#cd /hpcfs/home/ing_civil/da.martinez33/RENOBO/Proyectos/calculadora_financiera
#python download_SECOPII_table.py
cd /hpcfs/home/ing_civil/da.martinez33/UBC/6KAD5R5/SHE35318.20260421/20260417_LH00403_0172_A23MMVWLT4
module load fastqc
fastqc ./*.fastq.gz -o ./qualityControl/

echo -e "Finalicé la ejecución del script \n"
########################################################################################


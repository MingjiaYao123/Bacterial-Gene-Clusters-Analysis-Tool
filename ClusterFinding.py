
# Import parts of Biopython
from Bio import SeqIO
import difflib
import sys


# File path to your FASTA file
path_to_file=input("enter file name：")
theseq=input("enter the protein sequence of the core gene: ")
simi=int(input("set the good similarity for protein sequence comparing(30 is suggested)："))
cluster_length=int(input("set the length(bp) of each cluster(15000 is suggested): "))
output_file=open('gene_cluster.txt','w')
output_file.write('This output file contains the information of the found core genes and all the genes inside the found clusters from file '+path_to_file+',the size of each cluster is '+str(cluster_length)+' bp'+'\n')
n=0

# Open file with "with" statement to avoid problems with access 
# to original file (in case computer hangs
# or there will be any other problem)
with open(path_to_file, mode='r') as handle:

    # Use Biopython's parse function to process individual
    # FASTA records (thus reducing memory footprint)
    for record in SeqIO.parse(handle, 'fasta'):

        # Extract individual parts of the FASTA record
        identifier = record.id
        description = record.description
        sequence = record.seq
        
        lenth=len(theseq)
        exist=False
        s=sequence
        while len(s)>0:
            pos=min(s.find('ATG'),s.find('GTG'),s.find('TTG'))
            if pos==-1:
                break
            else:
                s=s[pos:]
                protseq=s.translate(to_stop=True)
                seq=difflib.SequenceMatcher(isjunk=None,a=theseq,b=protseq)
                diff=seq.ratio()*100
                diff=round(diff,1)
                if diff>=simi:
                    exist=True
                    break
                s=s[2:]
        
        if exist==True:
            n=n+1
            print('**************************************************************************************************')
            output_file.write('**************************************************************************************************'+'\n')
            print("Gene cluser "+ str(n)+", the core gene is:")
            output_file.write("Gene cluser "+ str(n)+", the core gene is:"+'\n')
            print('Processing the record {}:'.format(identifier))
            output_file.write('Processing the record {}:'.format(identifier)+'\n')
            print('Its description is: \n{}'.format(description))
            output_file.write('Its description is: \n{}'.format(description)+'\n')
            amount_of_nucleotides = len(sequence)
            print(sequence)
            output_file.write(str(sequence)+'\n')
            print(protseq)
            output_file.write(str(protseq)+'\n')
            print('Its sequence contains {} nucleotides.'.format(amount_of_nucleotides))
            output_file.write('Its sequence contains {} nucleotides.'.format(amount_of_nucleotides)+'\n')
            print("The similarity between this gene and input core gene is: "+str(diff)+'%')
            output_file.write("The similarity between this gene and input core gene is: "+str(diff)+'%'+'\n')
            cstart=0
            cend=0
            if not 'complement' in description and not 'join' in description and not '=<' in description and not '..>' in description:
                cpos1=description.find('location=')
                cloc1=cpos1+9
                cpos2=description.find('..')
                cloc2=cpos2+2
                cpos3=description.find('] [gbkey=CDS')
                cstart=int(description[cloc1:cpos2])
                cend=int(description[cloc2:cpos3])
            if 'complement' in description and not 'join' in description and not '=<' in description and not '..>' in description:
                cpos1=description.find('location=')
                cloc1=cpos1+20
                cpos2=description.find('..')
                cloc2=cpos2+2
                cpos3=description.find(')] [gbkey=CDS')
                cstart=int(description[cloc1:cpos2])
                cend=int(description[cloc2:cpos3])
            with open(path_to_file, mode='r') as handle2:
                n2=0
                for record2 in SeqIO.parse(handle2, 'fasta'):
                    identifier2 = record2.id
                    description2 = record2.description
                    sequence2 = record2.seq
                    
                    if not 'complement' in description2 and not 'join' in description2 and not '=<' in description2 and not '..>' in description2:
                        inside=False
                        pos1=description2.find('location=')
                        loc1=pos1+9
                        pos2=description2.find('..')
                        loc2=pos2+2
                        pos3=description2.find('] [gbkey=CDS')
                        start=int(description2[loc1:pos2])
                        end=int(description2[loc2:pos3])
                        if end >= (cend-cluster_length)and start <= (cstart+cluster_length):
                            inside=True
                        if inside==True:
                            n2=n2+1
                            print('----------------------------------------------------------')
                            output_file.write('----------------------------------------------------------'+'\n')
                            print('Cluster '+ str(n) +',gene '+str(n2))
                            output_file.write('Cluster '+ str(n) +',gene '+str(n2)+'\n')
                            print('Processing the record {}:'.format(identifier2))
                            output_file.write('Processing the record {}:'.format(identifier2)+'\n')
                            print('Its description is: \n{}'.format(description2))
                            output_file.write('Its description is: \n{}'.format(description2)+'\n')
                            amount_of_nucleotides2 = len(sequence2)
                            print(sequence2)
                            output_file.write(str(sequence2)+'\n')
                            print('Its sequence contains {} nucleotides.'.format(amount_of_nucleotides2))
                            output_file.write('Its sequence contains {} nucleotides.'.format(amount_of_nucleotides2)+'\n')

                    if 'complement' in description2 and not 'join' in description2 and not '=<' in description2 and not '..>' in description2:
                        inside=False
                        pos1=description2.find('location=')
                        loc1=pos1+20
                        pos2=description2.find('..')
                        loc2=pos2+2
                        pos3=description2.find(')] [gbkey=CDS')
                        start=int(description2[loc1:pos2])
                        end=int(description2[loc2:pos3])
                        if end >= (cend-cluster_length)and start <= (cstart+cluster_length):
                            inside=True
                        if inside==True:
                            n2=n2+1
                            print('----------------------------------------------------------')
                            output_file.write('----------------------------------------------------------'+'\n')
                            print('Cluster '+ str(n) +',gene '+str(n2))
                            output_file.write('Cluster '+ str(n) +',gene '+str(n2)+'\n')
                            print('Processing the record {}:'.format(identifier2))
                            output_file.write('Processing the record {}:'.format(identifier2)+'\n')
                            print('Its description is: \n{}'.format(description2))
                            output_file.write('Its description is: \n{}'.format(description2)+'\n')
                            amount_of_nucleotides2 = len(sequence2)
                            print(sequence2)
                            output_file.write(str(sequence2)+'\n')
                            print('Its sequence contains {} nucleotides.'.format(amount_of_nucleotides2))
                            output_file.write('Its sequence contains {} nucleotides.'.format(amount_of_nucleotides2)+'\n')
                print('totally '+str(n2)+' genes in this cluster')
                output_file.write('totally '+str(n2)+' genes in this cluster'+'\n')
                
output_file.close


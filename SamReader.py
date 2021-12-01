#!/usr/bin/python3
#-*- coding : utf-8 -*-


#__authors__ = ("Thibaud MARIN", "Kélian PEREZ")
#__contact__ = ("thibaud.marin@etu.umontpellier.fr","kelian.perez@etu.umontpellier.fr")
#__version__ = "1.0"
#__date__ = "12/14/2021"
#__licence__ ="This program is free software: you can redistribute it and/or modify
        #it under the terms of the GNU General Public License as published by
        #the Free Software Foundation, either version 3 of the License, or
        #(at your option) any later version.
        #This program is distributed in the hope that it will be useful,
        #but WITHOUT ANY WARRANTY; without even the implied warranty of
        #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        #GNU General Public License for more details.
        #You should have received a copy of the GNU General Public License
       # along with this program. If not, see <https://www.gnu.org/licenses/>."


     
    ### OPTION LIST:
        ##-h or --help : help information
        ##-i or --input: input file (.sam)
        ##-o or --output: output name files (.txt)

    #Synopsis:
        ##SamReader.py -h or --help # launch the help.
        ##SamReader.py -i or --input <file> # Launch SamReader to analyze a samtools file (.sam) and print the result in the terminal
        ##SamReader.py -i or --input <file> -o or --output <name> # Launch SamReader to analyze a samtools file (.sam) and print the result in the file called <name>
  


############### IMPORT MODULES ###############

import os, sys, re


############### FUNCTIONS TO :

## 1/ Check,
import os, sys
def testFile(fichier_entree) :
    file_is_correct=True
    if os.path.exists(fichier_entree):
        if os.path.isfile(fichier_entree):
            if os.stat(fichier_entree).st_size != 0:
                fichier = open(fichier_entree, "r")
                cpt=0
                for line in fichier:
        
                    if not line.startswith("@"):
                        cpt+=1
                        col_line=line.split('\t')
                        if len(col_line)>=11 :## Test du nbe de colonnes
                            if (int(col_line[1])>(2^12-1) or int(col_line[1])<0): ## Test du Flag
                                if re.match(r"\*|([0-9]+[MIDNSHPX=])+",col_line[5]):## Test du CIGAR
                                    pass
                                else:
                                    print("ERROR_FILE: The CIGAR value: "+col_line[5]+", is not valid.")
                                    file_is_correct=False
                                    break
                            else:
                                print("ERROR_FILE: The FLAG value, is out of the expected values [0,4095]")
                                file_is_correct=False
                                break
                        else :
                            print("ERROR_FILE: The number of columns is under the minimum required for a SAM file.")
                            file_is_correct=False
                            break
                        if (cpt>100):
                            break
            else:
                print('ERROR_FILE: '+fichier_entree+' is empty')
                file_is_correct=False
        else:
            print('ERROR_USER: Please insert a file')
            file_is_correct=False
    else:
        print('ERROR_PATH: Entry does not exists')
        file_is_correct=False
    if file_is_correct:
        print("Entry file is a correct SAM")

## 2/ Read, 



## 3/ Store,

def store_sam(file_input):
    ## Create a dictionnary with FLAG as keys and a list of important informations for each line (Name(col_line[0]), FLAG(col_line[1]), CIGAR(col_line[5]), sequence(col_line[9]))
    dico_sam={}
    
    ## Open the input file with "read" option
    with open(file_input,"r") as file_lines:
        for line in file_lines:
            if not line.startswith("@"):
                col_line = line.split("\t")
                flag = col_line[1]
                major_element_line=col_line[0]+"\t"+col_line[1]+"\t"+col_line[5]+"\t"+col_line[9]
                if flag in dico_sam:
                    dico_sam[flag].append(major_element_line) ## Attention, il faut choisir les element à écrire dans la liste du dico !!
                else:
                    dico_sam[flag]=[major_element_line]
    return dico_sam

## 4/ Analyse 










            
#### Convert the flag into binary ####
def flagBinary(flag) :

    flagB = bin(int(flag)) # Transform the integer into a binary.
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 
    if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
        add = 12 - len(flagB) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
        for t in range(add):
            flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
    return flagB















#### Analyze the unmapped reads (not paired) ####
def unmapped(dico_sam):
    
    unmapped_count = 0
    with open ("only_unmapped.fasta", "a+") as unmapped_fasta, open("summary_unmapped.txt", "w") as summary_file:
        for key_dico in dico_sam: # Parse the keys of the dictionnary, in our case all the possible Flag
            flag = flagBinary(key_dico) # Transform the flag into binary

            if int(flag[-3]) == 1: # Check if the condition 3 is True: 'Read unmapped'
                unmapped_count += len(dico_sam[key_dico])      # add the number of element in the list for the corresponding Flag (= line number with this flag)
                for line in dico_sam[key_dico]:                # Parse the lines having the corresponding flag
                    unmapped_fasta.write(str(line+"\n"))            # Write the line into a file called 'only_unmapped.fasta'

        summary_file.write("Total unmapped reads: " + str(unmapped_count) + "\n") # Write the total number of unmapped reads into a file called 'summary_unmapped.txt'
        return unmapped_count

    







    
#### Analyze the partially mapped reads ####
def partiallyMapped(dico_sam):
    
    partially_mapped_count = 0

    with open ("only_partially_mapped.fasta", "a+") as partially_mapped_fasta, open("summary_partially_mapped.txt", "w") as summary_file:
        for key_dico in dico_sam:
            flag = flagBinary(key_dico) # We compute the same
            for line in dico_sam[key_dico]:
                col_line=line.split('\t')
                cigar_list_dico=col_line[2]
                if int(flag[-2]) == 1: 
                    if cigar_list_dico != "100M":
                        partially_mapped_count += 1
                        partially_mapped_fasta.write(str(line+"\n"))

        summary_file.write("Total partially mapped reads: " + str(partially_mapped_count) + "\n") 
        return partially_mapped_count

    
#### Analyze pair of reads where read is mapped and mate unmapped #### -> check du FLAG pour mate unmapped et CIGAR = 100M
def Mapped_Unmapped(dico_sam):
    read_mapped_mate_unmapped_count = 0

    with open ("read_mapped_mate_unmapped.fasta", "a+") as read_mapped_mate_unmapped_fasta, open("summary_read_mapped_mate_unmapped.txt", "w") as summary_file:
        for key_dico in dico_sam:
            flag = flagBinary(key_dico) # We compute the same
            for line in dico_sam[key_dico]:
                col_line=line.split('\t')
                cigar_list_dico=col_line[2]
                if int(flag[-4]) ==1:
                    if re.match(r"\*|([0-9]+[M=])+",cigar_list_dico) :
                        read_mapped_mate_unmapped_count += 1
                        read_mapped_mate_unmapped_fasta.write(str(line+"\n"))

        summary_file.write("Total pair of reads where a read is mapped and his mate is unmapped: " + str(read_mapped_mate_unmapped_count) + "\n") 
        return read_mapped_mate_unmapped_count

#### Fonction that return the mate flag ####

def flag_mate(flag):
    flag_bin=flagBinary(flag)
    flag_mate_bin=flag_bin
    flag_mate=0
    if flag_bin[-5]=='1': # Modify the sense
        flag_mate_bin[-5]='0'
        flag_mate_bin[-6]='1'
    else:
        flag_mate_bin[-5]='1'
        flag_mate_bin[-6]='0'

    if flag_mate_bin[-7]=='1': # Modify the first/secon in pair
        flag_mate_bin[-7]='0'
        flag_mate_bin[-8]='1'
    else:
        flag_mate_bin[-7]='1'
        flag_mate_bin[-8]='0'
        
    j=0
    for i in range(len(flag_mate_bin),0,-1):
        if flag_mate_bin[i-1]=='1':
            flag_mate+=2**j
        j+=1
    
    return flag_mate



#### Analyze the reads where one is mapped and the mate is partially mapped (Using Flag and CIGAR)####

def one_partially_mapped(dico_sam):
    pair_one_partially_mapped_count = 0
    
    with open ("one_partially_mapped.fasta", "a+") as one_partially_mapped_fasta, open("summary_one_partially_mapped.txt", "w") as summary_file:
        for key_dico in dico_sam:
            flag = flagBinary(key_dico) # We compute the same
            if flag[-2]=='1':                                   # Check if the 2nd element in the FLAG is 1 (paired in proper pair)
                for line in dico_sam[key_dico]:                 # Parse the lines with the corresponding flag
                    col_line=line.split('\t')                      
                    cigar=col_line[2]                           # CIGAR is stocked into the variable cigar
                    if not re.match(r"\*|([0-9]+[M=])+",cigar): # Check if the read is partially mapped
                        name_read=col_line[0]
                        
                       
                        for line_mate in dico_sam[str(flag_mate(key_dico))]:                        # Parse lines with complementary flag, to find the mate
                            col_line_mate=line_mate.split('\t')
                            cigar_mate=col_line_mate[2]


                            if ((col_line_mate[0]==name_read) & (bool(re.match(r"\*|([0-9]+[M=])+",cigar_mate)))):  # Check: if same name as mate, if its CIGAR is full M (full mapped)
                                pair_one_partially_mapped_count += 1                                                
                                one_partially_mapped_fasta.write(str(line+"\n"))                                    # Write the couple of reads in the file
                                one_partially_mapped_fasta.write(str(line_mate+"\n"))
                                break
                        
                        

        summary_file.write("Total pair of reads with one partially mapped and one mapped correctly: " + str(pair_one_partially_mapped_count)+ "\n") 
        return pair_one_partially_mapped_count

    

    
### Analyse the CIGAR = regular expression that summarise each read alignment ###
def readCigar(cigar): 
   
    ext = re.findall('\w',cigar) # split cigar 
    key=[] 
    value=[]    
    val=""

    for i in range(0,len(ext)): # For each numeric values or alpha numeric
        if (ext[i] == 'M' or ext[i] == 'I' or ext[i] == 'D' or ext[i] == 'S' or ext[i] == 'H' or ext[i] == "N" or ext[i] == 'P'or ext[i] == 'X'or ext[i] == '=') :
            key.append(ext[i])
            value.append(val)
            val = ""
        else :
            val = "" + val + ext[i]  # Else concatenate in order of arrival
    
    dico = {}
    n = 0
    for k in key:   # Dictionnary contruction in range size lists              
        if k not in dico.keys():    # for each key, insert int value
            dico[k] = int(value[n])   # if key not exist, create and add value
            n += 1
        else:
            dico[k] += int(value[n])  # inf key exist add value
            n += 1
    return dico

### Analyse the CIGAR = regular expression that summarise each read alignment ###
def percentMutation(dico):
        
    totalValue = 0 # Total number of mutations
    for v in dico :
        totalValue += dico[v]

    mutList = ['M','I','D','S','H','N','P','X','=']
    res = ""
    for mut in mutList : # Calculated percent of mutation if mut present in the dictionnary, else, percent of mut = 0
        if mut in dico.keys() :
            res += (str(round((dico[mut] * 100) / totalValue, 2)) + ";")
        else :
            res += ("0.00" + ";")
    return res

def globalPercentCigar():
    """
      Global representation of cigar distribution.
    """
    
    with open ("outpuTable_cigar.txt","r") as outpuTable, open("Final_Cigar_table.txt", "w") as FinalCigar:
        nbReads, M, I, D, S, H, N, P, X, Egal = [0 for n in range(10)]

        for line in outpuTable :
            mutValues = line.split(";")
            nbReads += 2
            M += float(mutValues[2])+float(mutValues[12])
            I += float(mutValues[3])+float(mutValues[13])
            D += float(mutValues[4])+float(mutValues[14])
            S += float(mutValues[5])+float(mutValues[15])
            H += float(mutValues[6])+float(mutValues[16])
            N += float(mutValues[7])+float(mutValues[17])
            P += float(mutValues[8])+float(mutValues[18])
            X += float(mutValues[9])+float(mutValues[19])
            Egal += float(mutValues[10])+float(mutValues[20])

        FinalCigar.write("Global cigar mutation observed :"+"\n"
                        +"Alignlent Match : "+str(round(M/nbReads,2))+"\n"
                        +"Insertion : "+str(round(I/nbReads,2))+"\n"
                        +"Deletion : "+str(round(D/nbReads,2))+"\n"
                        +"Skipped region : "+str(round(S/nbReads,2))+"\n"
                        +"Soft Clipping : "+str(round(H/nbReads,2))+"\n"
                        +"Hard Clipping : "+str(round(N/nbReads,2))+"\n"
                        +"Padding : "+str(round(P/nbReads,2))+"\n"
                        +"Sequence Match : "+str(round(Egal/nbReads,2))+"\n"
                        +"Sequence Mismatch : "+str(round(X/nbReads,2))+"\n")


 
#### Summarise the results ####

#def Summary(fileName):
    
   
#### Help function ####

def help():
    print("\nAUTHORS: Thibaud MARIN, Kélian PEREZ\n\n"+
          "USAGE:\tSamReader.py [options] <file_in.sam>|<file_out.txt>\n\n"+
          "OPTIONS:\n -h|--help\t  \t Display AUTHORS, USAGE and OPTIONS for SamReader.py\n"+
          " -i|--input\tFILE\t Input file path\n"+
          " -o|--output\tFILE\t Output file name [file_out.txt]\n\n")

#### Main function ####

def main(argv):
    #fichier_entree=sys.argv[1]
    #testFile(fichier_entree)
    dico_files={
        "file_in":"file_in",
        "file_out":"./file_out.txt",
    }
    
    for i in range(0,len(sys.argv),1):
        if ( sys.argv[i] == "-i" ) or ( sys.argv[i] == "--input" ):
            dico_files["file_in"]=sys.argv[i+1]
        elif ( sys.argv[i] == "-o" ) or ( sys.argv[i] == "--output" ):
            dico_files["file_out"]=sys.argv[i+1]
        elif ( sys.argv[i] == "-h" ) or ( sys.argv[i] == "--help" ) or ( len(sys.argv) == 1 ):
            help()

    dico_file_sam=store_sam(dico_files["file_in"])
    #unmapped(dico_file_sam)
    #partiallyMapped(dico_file_sam)
    #testFile(dico_files["file_in"])
    #checkUtf8((dico_files["file_in"])
    Mapped_Unmapped(dico_file_sam)
    
############### LAUNCH THE SCRIPT ###############

if __name__ == "__main__":
    main(sys.argv[1:])

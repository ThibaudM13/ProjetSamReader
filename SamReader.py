#!/usr/bin/python3
#-*- coding : utf-8 -*-


#__authors__ = ("Thibaud MARIN", "Kélian PEREZ")
#__contact__ = ("thibaud.marin@etu.umontpellier.fr","kelian.perez@etu.umontpellier.fr")
#__version__ = "1.2"
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
def testFile(given_file) :
    print("Testing file: ",end='')

    file_is_correct=True
    nb_header_lines=0
    
    ## Tests on given file
    
    if os.path.exists(given_file):
        if os.path.isfile(given_file):
            if os.stat(given_file).st_size != 0:
                if (((given_file.split('.'))[-1])=="sam"):
                    fichier = open(given_file, "r")
                    cpt=0
                    for line in fichier:
                        if line.startswith("@"): ## count the header lines
                            nb_header_lines+=1
                        else:
                            cpt+=1
                            col_line=line.split('\t')
                            
                            ## Test du nbe de colonnes
                            if len(col_line)>=11 :
                                if (int(col_line[1])>(2^12-1) or int(col_line[1])<0): ## Flag test
                                    if re.match(r"\*|([0-9]+[MIDNSHPX=])+",col_line[5]):## CIGAR test
                                        pass
                                    else:
                                        print(f"ERROR_FILE: The CIGAR value: {col_line[5]}, is not valid.")
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
                            if (cpt=100): # Check the format and the values of each field 
                                break
                else:
                    print("ERROR_USER: extension different from '.sam', please try again with a SAM file.")            
                    file_is_correct=False
            else:
                print(f"ERROR_FILE: {given_file} is empty")
                file_is_correct=False
        else:
            print('ERROR_USER: Please insert a file')
            file_is_correct=False
    else:
        print('ERROR_PATH: Entry does not exists')
        file_is_correct=False
        


    ## Bilan
    print("done.")
    if ((cpt==100) and file_is_correct):
        return True
    else:
        return False

        
## 2/ Read, 



## 3/ Store,

def storeSam(file_input):
    print("Storing file into dictionnary: ",end='')
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
    print("done.")
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
def unmapped(dico_sam, file_out):
    print("Function unmapped: ",end='')
    unmapped_count = 0
    with open ("only_unmapped.fasta", "w") as unmapped_fasta, open(file_out, "a+") as summary_file:
        for key_dico in dico_sam: # Parse the keys of the dictionnary, in our case all the possible Flag
            flag = flagBinary(key_dico) # Transform the flag into binary

            if int(flag[-3]) == 1: # Check if the condition 3 is True: 'Read unmapped'
                unmapped_count += len(dico_sam[key_dico])      # add the number of element in the list for the corresponding Flag (= line number with this flag)
                for line in dico_sam[key_dico]:                # Parse the lines having the corresponding flag
                    col_line=line.split('\t')
                    unmapped_fasta.write(f">{col_line[0]} function:unmapped\n{col_line[3]}\n")           # Write the line into a file called 'only_unmapped.fasta'

        summary_file.write(f"\nTotal unmapped reads: {unmapped_count}\n") # Write the total number of unmapped reads into a file summary.
    print("done.")

    




    
#### Analyze the partially mapped reads ####
def partiallyMapped(dico_sam, file_out):
    print("Function partiallyMapped: ",end='')
    partially_mapped_count = 0

    with open ("partially_mapped.fasta", "w") as partially_mapped_fasta, open(file_out, "a+") as summary_file:
        for key_dico in dico_sam:
            flag = flagBinary(key_dico)
            for line in dico_sam[key_dico]:
                col_line=line.split('\t')
                cigar=col_line[2]
                if re.fullmatch(r".*\d[SH].*",cigar): # Selection of all reads partially mapped (Cigar different of only M)
                    partially_mapped_count += 1
                    partially_mapped_fasta.write(f">{col_line[0]} function:partiallyMapped\n{col_line[3]}\n")

        summary_file.write(f"Total partially mapped reads: {partially_mapped_count}\n") 
    print("done.")




    
#### Analyze pair of reads where read is mapped and mate unmapped #### -> check du FLAG pour mate unmapped et CIGAR = 100M
def mappedUnmapped(dico_sam, file_out):
    print("Function Mapped_Unmapped: ",end='')
    read_mapped_mate_unmapped_count = 0

    with open ("read_mapped_mate_unmapped.fasta", "w") as read_mapped_mate_unmapped_fasta, open(file_out, "a+") as summary_file:
        for key_dico in dico_sam:
            flag = flagBinary(key_dico)
            if (flag[-4] =='1' and flag[-3] =='0'):
                for line in dico_sam[key_dico]:
                    col_line=line.split('\t')
                    cigar=col_line[2]
                    if re.fullmatch(r"\d*[M]",cigar):
                        name_read=col_line[0]
                        for line_mate in dico_sam[str(flagMate(key_dico))]:                        # Parse lines with complementary flag, to find the mate
                            col_line_mate=line_mate.split('\t')
                            cigar_mate=col_line_mate[2]
                            if (col_line_mate[0]==name_read):                                       # Check: if same name as mate
                                read_mapped_mate_unmapped_count += 1
                                read_mapped_mate_unmapped_fasta.write(f">{col_line[0]} function:Mapped_Unmapped read\n{col_line[3]}\n")               # Write the couple of reads in the file
                                read_mapped_mate_unmapped_fasta.write(f">{col_line_mate[0]} function:Mapped_Unmapped mate\n{col_line_mate[3]}\n")
                                break

        summary_file.write(f"Total pair of reads where a read is mapped and his mate is unmapped: {read_mapped_mate_unmapped_count}\n") 
    print("done.")
    
   
   
   
   
#### Fonction that return the mate flag ####

def flagMate(flag):
    flag_bin=flagBinary(flag)
    flag_Mate_bin=flag_bin
    flag_Mate=0
    if flag_bin[-3]=='1' and flag_bin[-4]=='0': # Modify the read/ mate unmapped
        flag_Mate_bin[-3]='0'
        flag_Mate_bin[-4]='1'
    elif flag_bin[-4]=='1'and flag_bin[-3]=='0':
        flag_Mate_bin[-3]='1'
        flag_Mate_bin[-4]='0'
 
    if flag_bin[-5]=='1' and flag_bin[-6]=='0': # Modify the sense
        flag_Mate_bin[-5]='0'
        flag_Mate_bin[-6]='1'
    elif flag_bin[-6]=='1' and flag_bin[-5]=='0':
        flag_Mate_bin[-5]='1' 
        flag_Mate_bin[-6]='0'

    if flag_bin[-7]=='1' and flag_bin[-8]=='0': # Modify the first/second in pair
        flag_Mate_bin[-7]='0'
        flag_Mate_bin[-8]='1'
    elif flag_bin[-8]=='1' and flag_bin[-7]=='0':
        flag_Mate_bin[-7]='1'
        flag_Mate_bin[-8]='0'
        
    j=0
    for i in range(len(flag_Mate_bin),0,-1):
        if flag_Mate_bin[i-1]=='1':
            flag_Mate+=2**j
        j+=1
    
    return flag_Mate



#### Analyze the reads where one is mapped and the mate is partially mapped (Using Flag and CIGAR)####

def onePartiallyMapped(dico_sam,file_out):
    print("Function one_partially_mapped: ",end='')
    pair_one_partially_mapped_count = 0
    
    with open ("one_partially_mapped.fasta", "w") as one_partially_mapped_fasta, open(file_out, "a+") as summary_file:
        for key_dico in dico_sam:
            flag = flagBinary(key_dico) # We compute the same
            if flag[-2]=='1':                                   # Check if the 2nd element in the FLAG is 1 (paired in proper pair)
                for line in dico_sam[key_dico]:                 # Parse the lines with the corresponding flag
                    col_line=line.split('\t')                      
                    cigar=col_line[2]                           # CIGAR is stocked into the variable cigar
                    if re.fullmatch(r".*\d[SH].*",cigar): # Check if the read is partially mapped
                    
                        for line_mate in dico_sam[str(flagMate(key_dico))]:                        # Parse lines with complementary flag, to find the mate
                            col_line_mate=line_mate.split('\t')
                            cigar_mate=col_line_mate[2]
                            
                            if ((col_line_mate[0]==col_line[0]) & (bool(re.fullmatch(r"\d*[M]",cigar_mate)))):  # Check: if same name as mate, if its CIGAR is full M (full mapped)
                                pair_one_partially_mapped_count += 1                                                
                                one_partially_mapped_fasta.write(f">{col_line[0]} function:one_partially_mapped read partially mapped\n{col_line[3]}\n")       # Write the couple of reads in the file
                                one_partially_mapped_fasta.write(f">{col_line_mate[0]} function:one_partially_mapped mate mapped\n{col_line_mate[3]}\n")
                                break
                        
                        

        summary_file.write(f"Total pair of reads with one partially mapped and one mapped correctly: {pair_one_partially_mapped_count}\n") 
    print("done.")

    

    
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



def globalPercentCigar(dico_sam,file_out):
    """
        Preparation of the outpuTable_cigar.txt file.
    """
    print("Function globalPercentCigar (preparation): ",end='')
    with open ("outpuTable_cigar.txt", "w") as outputTable:
        for flag in dico_sam:
            for line in dico_sam[flag]:
                col_line=line.split('\t')
                outputTable.write(f"{percentMutation(readCigar(col_line[2]))}\n")


    print("done.")
    """
      Global representation of cigar distribution.
    """
    print("Function globalPercentCigar (analyse and writing): ",end='')
    with open ("outpuTable_cigar.txt","r") as outpuTable, open(file_out, "a+") as FinalCigar:
        nbReads, M, I, D, S, H, N, P, X, Egal = [0 for n in range(10)]

        for line in outpuTable :
            mutValues = line.split(";")
            nbReads += 1
            M += float(mutValues[0])
            I += float(mutValues[1])
            D += float(mutValues[2])
            S += float(mutValues[3])
            H += float(mutValues[4])
            N += float(mutValues[5])
            P += float(mutValues[6])
            X += float(mutValues[7])
            Egal += float(mutValues[8])

        FinalCigar.write("\n\n\n\tGlobal cigar mutation observed :"+"\n\n"
                        +"Alignement Match : "+str(round(M/nbReads,2))+"\n"
                        +"Insertion : "+str(round(I/nbReads,2))+"\n"
                        +"Deletion : "+str(round(D/nbReads,2))+"\n"
                        +"Skipped region : "+str(round(S/nbReads,2))+"\n"
                        +"Soft Clipping : "+str(round(H/nbReads,2))+"\n"
                        +"Hard Clipping : "+str(round(N/nbReads,2))+"\n"
                        +"Padding : "+str(round(P/nbReads,2))+"\n"
                        +"Sequence Match : "+str(round(Egal/nbReads,2))+"\n"
                        +"Sequence Mismatch : "+str(round(X/nbReads,2))+"\n")
    os.remove("outpuTable_cigar.txt")
    print("done.")




 
#### Summarise the results ####

def summary(dico_sam, fileSummaryName):
    with open(fileSummaryName,"w") as f:
        f.write(f"\t\t\t\t   GLOBAL SUMMARY OF YOUR SAM FILE \n"+
                "\tNumber of reads with characteristics:\n")
        
    unmapped(dico_sam, fileSummaryName)
    mappedUnmapped(dico_sam, fileSummaryName)
        
    partiallyMapped(dico_sam,fileSummaryName)
    onePartiallyMapped(dico_sam,fileSummaryName)
    
    globalPercentCigar(dico_sam,fileSummaryName)
   

   
#### Help function ####

def help():
    print("\nAUTHORS: Thibaud MARIN, Kélian PEREZ\n\n"+
          "USAGE:\tSamReader.py [options] <file_in.sam>|<file_out.txt>\n\n"+
          "OPTIONS:\n -h|--help\t  \t Display AUTHORS, USAGE and OPTIONS for SamReader.py\n"+
          " -i|--input\tFILE\t Input file path\n"+
          " -o|--output\tFILE\t Output file name [file_out.txt]\n\n")




#### Main function ####

def main(argv):
    dico_files={
        "file_in":"file_in",
        "file_out":"./file_out_summary.txt",
    }

    
    ## Arguments gestion
    for i in range(0,len(sys.argv),1):
        liste_arg=["-i","--input","-o","--output","-h","--help"]
        
        if ( sys.argv[i] in liste_arg[0:2:1] ):
        
            # Check if a file name/path is given after the option -i|--input and different from an option
            if i+1<len(sys.argv) and sys.argv[i+1] not in liste_arg:
                dico_files["file_in"]=sys.argv[i+1]
            else:
                print("Please insert a file name/path after the -i|--input option")
                exit()
        
        elif ( sys.argv[i] in liste_arg[2::1][:2:1] ):
        
            # Check if a file name or path is given after the option -o|--output and different from an option
            if i+1<len(sys.argv) and sys.argv[i+1] not in liste_arg :
                dico_files["file_out"]=sys.argv[i+1]
            else:
                print("Please insert a file name/path after the -o|--output option")
                exit()
        
        elif ( sys.argv[i] in liste_arg[4:2:1][:2:1] ) or ( len(sys.argv) == 1 ):
            help()
            exit()
        
        elif ((sys.argv[i] not in liste_arg) and i!=0 and (sys.argv[i-1] not in liste_arg)):
            print(f"Please enter correct arguments (for more information, use the --help|-h option).")
            exit()


    ## Launch the script Summary, with all fonctions if the file is correct
    if testFile(dico_files["file_in"]):
        print("Entry file is a correct SAM")
        dico_sam=storeSam(dico_files["file_in"])
        summary(dico_sam,dico_files["file_out"])
    
    
############### LAUNCH THE SCRIPT ###############

if __name__ == "__main__":
    main(sys.argv[1:])

<div id="top"></div>


<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/ThibaudM13/ProjetSamReader">
    <img src="Logo_samreader.png">
  </a>

<h3 align="center">Projet SamReader</h3>

  <p align="center">
    Parse and give basic informations about a SAM file.
    <br />
    <a href="https://github.com/ThibaudM13/ProjetSamReader"><strong>Explore the GitHub »</strong></a>
    <br />
    <br />
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

In order to run optimally our SamReader.py, this are configuration and packages our script uses:

* Python 3
  ```sh
  # If your Python version is under 3.0, the script could not work, please install Python 3.0:
  
  sudo apt install python3
  ```
  
* Packages used: re, os, sys (basic Python packages)

### Installation

Download the script using one of the following way:

  ```sh
  # By using wget
  wget https://github.com/ThibaudM13/ProjetSamReader/blob/master/SamReader.py
  
  
  # By using git clone (all elements will be downloaded)
  git clone https://github.com/ThibaudM13/ProjetSamReader.git
  ```


<!-- USAGE EXAMPLES -->
## Usage

### Input

SAM file with '.sam' extension (Be sure that the script will verify the input format)

### Launch the script

```sh
./SamReader.py  --input|-i /your/path/to/your/SAM/file.sam 
                < --output|-o /your/path/to/your/summary_file.txt >
                < --help|-h >
```

### Output

A summary file (by default named 'file_out_summary.txt', on your current repertory) containing:

* Number of reads with characteristics targeted by the functions: <br />
  - reads unmapped
  - pair where a read is correctly mapped (CIGAR= 100M), and the mate is unmapped
  - reads partially mapped (CIGAR with soft clipped)
  - pair when a read is correctly mapped (CIGAR= 100M), and the mate is partially mapped
  <br/>
  <br/>

* Global characteristics about the mapping, based on the CIGAR

<br/>
  
Fasta files (one file for each function)

- File name:

  `<function_name>.fasta`

- File content:

  * ID and comments <br />
    format: `> <read_name> function:<function_name> <read_characteristic>`
  
  * Sequence 




<!-- LICENSE -->
## License

This program is free software: you can **redistribute** it and/or **modify** it under the terms of the **GNU General Public License** as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

See the GNU General Public License for more details ([GNU website](https://www.gnu.org/licenses/)).



<!-- CONTACT -->
## Contact

MARIN Thibaud - [thibaud.marin@etu.umontpellier.fr](mailto:thibaud.marin@etu.umontpellier.fr) - [Linkedin](https://www.linkedin.com/in/thibaud-marin/)

PEREZ Kélian - [kelian.perez@etu.umontpellier.fr](mailto:kelian.perez@etu.umontpellier.fr) - [Linkedin](https://www.linkedin.com/in/k%C3%A9lian-perez-b0217a208/)

<p align="right">(<a href="#top">back to top</a>)</p>

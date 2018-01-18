# Linux Basics
Linux has been used for a large-scale collaborated research work environment. Many bioinformatics software has developed in Linux. In this section, we would like to get familiar with Linux operating system that requires minimal knowledge to run bioinformatics software tools properly.

## Setup and Requirements
### 1. Account Information
1. This is our first Bioinformatics workshop.
1. We will cover many basics of Linux to run bioinformatics software
1. No prerequisite is required.
1. For those of you registered for this workshop, we already created a Linux account and in advance sent to each of you credential information to access a Linux network at LRI (Lerner Research Institute).
    1. Linux login account ID/Password
    1. The Linux hostname you want to connect
        1. e.g., `scotch.lerner.ccf.org`

### 2. Hardwares
1. Laptop, Desktop, or any device (ssh terminal is available) that Cleveland Clinic approves to connect CCF network
1. Make sure that your computer has CCF WI-FI connection (not connect to *guest SSID)

### SSH terminal and X11 forwarding
#### 3. Only Windows user
1. X11 forwarding
    1. Download X11 forward client, [Xming](https://sourceforge.net/projects/xming/files/Xming/6.9.0.31/Xming-6-9-0-31-setup.exe/download)
    1. Install the software into your local directory (for example, `C:\Users\<username>\Downloads\Xming`)
    1. Open the file and click all ‘Yes’ or ‘OK’ to complete the installation
1. putty
    1. Download [putty](https://the.earth.li/~sgtatham/putty/latest/w32/putty.exe) into `C:\Users\<username>\Downloads\Xming`, the same directory you installed Xming above

### 4. Connect to LRI Linux node
1. 4 Linux servers are available but we will inform each of you only one IP address
    ```bash
    10.66.64.37
    10.66.64.38
    10.66.64.39
    10.66.64.40
    ```
1. [putty](images/windows_putty_login.png) for Windows users
    1. run putty
    1. type the IP address that we assign you above)        
    1. save a profile with a nickname (e.g., scotch) for the hostname
    1. click "open"
    1. type account ID (user name)
    1. type password
1. [teriminal](images/mac_os_x_terminal_login.png) for Mac OS users
    1. open `terminal`
    1. type `ssh -X username@<ip_address>` # <ip_address> that we assign you above
    1. type password

### 5. Code of conduct
1. Respect other users, don’t share your login account with the other, how to post a question/answer, and etc.
1. Refer [here](https://wiki.archlinux.org/index.php/Code_of_conduct) (3 min)

## Linux commands
### 6. The quickest way to learn Linux
1. Find an old laptop you don't use
1. Install [a Linux](https://distrowatch.com/) distro and play with it.

### 7. Help
1. You don't have to remember the Linux command all details
    1. The more you use, the quicker you get familiar with
    1. Help from Linux communities.
    1. There is a general pattern which appears in every command line. 
1. Interaction
    1. Try with simple inputs
    1. See message from terminal
    1. Try with more input options or complicated ones
    1. Follow an error message
1. Help
    1. type `man top` and enter
    1. `top -h` or `top --help`
    1. google search

### 8. Knowing about the Linux box you log in
1. Task: to figure out the Linux distribution
    1. `cat /proc/version`
    1. Cent OS, RedHat (rpm): Linux cluster environment
    1. Ubuntu, Debian (deb): Personal desktop/laptop
    1. SuSE Linux
1. Task: to figure out Linux kernel version
    1. `uname -a`
1. Task: to figure out the computer specification you currently access
    1. `cat /proc/cpuinfo`
    1. Q1: Figure out memory size
1. Task: to open process viewers and how many jobs are running
    1. `top`
    1. `htop`
1. Task: to know which file systems mounted and how much HDD spaces are available
    1. `df -h`
1. Q2: Can you repeat checking CPU information you typed in the task above?
    1. `ctrl+r` and type `info`
    1. `cat` and `pageup` or `pagedown`

### 9. Navigate directories and handling a file 
1. changing directory(`cd`), making directory (`mkdir`), and listing files (`ls`)
    1. `cd ~/` # Q2a: What is "~/"?
    1. `ls ~/`
    1. `ls -l ~/`
    1. `ls -la ~/` # How is it different from the other?
    1. `pwd` # What is the output?
    1. `mkdir tmp`
    1. `cd tmp`
    1. `ls -la ./` # Q3: What does `./` mean?
    1. `cd ../` # Q3a: What does `../` mean?
    
1. Copy and Paste
    1. `cd tmp`
    1. Q4: Figure out the current directory in a full pathname
    1. Double click the output to select the directory name
    1. For window user,
        1. type `ls -la`, space, and click right mouse button to paste the directory name  
    1. For Mac user,
        1. type `ls -la`, space, and click middle mouse button
1. File copy and view a file
    1. `cd ~/`
    1. `cp .bashrc ~/tmp/`
    1. Q5: change directory to `tmp`
    1. Q5a: view the file content of `.bashrc` using `cat` command
1. Delete the file
    1. `rm ~/tmp/.bashrc`

### 10. Environment Variables
1. Task: meet first default variable $HOME
    1. `ls -la ~/tmp` vs. `ls -la $HOME/tmp` #Q6: Is it same?
    1. `echo $HOME` # Append `$` at the beginning to access the variable content
1. Task: Define a variable
    1. `printenv` #List out all environment variables currently registered
    1. Q7: Assign `$HOME/tmp` with a new alias (e.g., TMP2)  
    1. Q8: Change the working directory to `$HOME/tmp` using the variable above
1. Task: Define an alias
    1. `alias l='ls --color=auto -lat'`
    1. `l ~/`

### 11. Script (File editor vi)
1. Task: Create your first shell script
    1. `cd ~/tmp`
    1. `vi hello_world.sh` #Q8a: is there any meaning `.sh` at the end of file? 
        1. type `i` #can you see "-- INSERT --" at the bottom left window?
        1. type `echo "Welecome Bioinformatics Workshop 2018, ${USER}"` at the first line
        1. type `echo "browsing ${HOME}/tmp"` in the next line
        1. type `ls -la ${HOME}/tmp` in the next line
    1. press `esc` key #"-- INSERT --" at the bottom left window should disappear
    1. `:`, `wq!`, and then hit `enter` #save and exit
    1. `cat hello_world.sh`
1. Task: Run the script
    1. `bash hello_world.sh`
1. Q8a: Play `vi` more
    1. `vi test.txt`
    1. Figure out how we can undo, redo, delete one line, search a word, and etc (3 min)
    1. refer to [vi manual](https://www.cs.umd.edu/~yhchan/vim.pdf) for more details

### 12. Prepare a project directory
1. Task: Create a project directory under your $HOME directory
    1. Q9: Change to `$HOME` directory
    1. Q10: Create a directory called `proojects` under the `$HOME` directory (I know there is a typo)
    1. Q11: Delete the typo directory
    1. Q12: Is it possible to restore the file or directory you mistakenly deleted?
    1. Q12a: Create a directory called `projects` under the `$HOME` directory
1. Task: Preparing the workshop sample files
    1. Q13: Copy our course material data file into `projects`
        1. `/Informatics_workshop/binf20180119_data.tar.gz` #TODO
        1. Use <tab> key to complete matched words
   1. Extract the tar.gz file
        1. `cd ~/projects`
        1. `tar zxvf ./bioinfo_2018.tar.gz`
   1. Task: Check the overall dirctory size
        1. type `du -hs ./bioinfo_2018`
   1. Q14: Backup. Create a new directory (e.g., `archives`) under `projects` and move the tar.gz file into the directory  

### 13. Search files
1. `cd ~/projects`
1. Scenario: a part of file name
    1. The file name contains "chr22" and file extension is ".gz". The file must be somewhere under `projects` 
    1. `find ./bioinfo_2018 -name "chr7*.gz" -type f` # man find
1. Scenario: a file containing a certain word
    1. The file extension is `.txt` and I remember that it contains a word like, `intergenic`.
    1. `grep -ilR --include=*.txt -e "intergenic" ./bioinfo_2018` # man grep

### 14. View a file and editing
1. Task: Count the total number of FASTQ files
    1. `cd bioinfo_2018`
    1. Q14a: Create an alias (e.g., `binf`) to change our project directory `~/projects/bioinfo_2018` (3 min)
    1. `find fastq -maxdepth 1 -type f | wc -l`
    1. `find fastq/*.fastq.gz -maxdepth 1 -type f | wc -l` 
        1. Only count `fastq.gz` files under `fastq` but not recursively
        1. Look for only file type not directory
        1. Q15: What is '`|`'?
        1. Q16: What is `wc -l`?
1. Task: Know who create a file and when it was modified
    1. `binf`
    1. `cd fastq`
    1. type `stat chr7.R1_001.fastq.gz`

1. Task: View a long wide TAB delimited table file
    1. `binf`
    1. `cd examples`
    1. `head -n 100 annovar.chr22.snps.hg19_multianno.txt`
    1. `tail -n 100 annovar.chr22.snps.hg19_multianno.txt`
    1. `less -S ./annovar.chr22.snps.hg19_multianno.txt` # to view the file and scroll up/down
    1. type `G` #to go the end of page
    1. type `/` # can you see `/` at the left bottom of window?
    1. type `PRAME`, press `n`, and press `N`
    1. type `g` #to go to the first page
    1. type `q` #to exit

### 15. Compress file or directory
1. a single file
    1. `gzip annovar.chr22.snps.hg19_multianno.txt`
    1. `less -S annovar.chr22.snps.hg19_multianno.txt.gz` # Can you view the compressed file content with `less`?
    1. `gunzip annovar.chr22.snps.hg19_multianno.txt`
1. archiving a directory    
    1. Q16a: Change to the parent directory
    1. `tar czvf examples.tar.gz ./examples`
    1. `rm -rf examples.tar.gz`

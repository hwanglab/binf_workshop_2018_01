#1. Customize your `bash` environment
1. Task: To configure your `bash` shell configuration to make life easier
    1. `cd ~/`
    1. `less ~/.bashrc`
    1. Q16b: Backup `.bashrc` using Linux copy command
    1. `vi .bashrc`
    1. add the following two lines at the end of file
    1. `alias l='ls --color=auto -lat'`
    1. `alias ll='ls --color=auto -lht'`
    1. `alias binf='cd $HOME/projects/bioinfo_2018'`
    1. `alias fastq='cd $HOME/projects/bioinfo_2018/fastq'`
    1. Q16c: Save the file and exit
1. Task: Configure `~/.bash_profile`
    1. Q16d: Make a directory, named `bin` under your $HOME directory
    1. `echo $PATH` to see if $HOME/bin is included in $PATH
    1. otherwise, open `~/.bash_profile` add the directory in $PATH
1. Task: Use Linux softlink commands, copy all bioinformatics program file into a directory in $PATH
    1. `binf`
    1. cd `utils`
    1. check if `lnS_path_setup.sh` exists
    1. `bash lnS_path_setup.sh`
```bash
#!/bin/bash -l

if [ ! -d ${HOME}/bin ]; then
    echo "Creating ${HOME}/bin directory ..."
    mkdir ${HOME}/bin
    echo "Done."
fi

ln -s /usr/local/tools/gatk/GenomeAnalysisTK.jar $HOME/bin/
ln -s /usr/local/tools/picard/picard.jar $HOME/bin/
ln -s /usr/local/tools/igv/igv.sh $HOME/bin/
ln -s /usr/local/tools/cutadapt/cutadapt $HOME/bin/
ln -s /usr/local/tools/bwa/bwa $HOME/bin/
ln -s /usr/local/tools/bedtools/bedtools $HOME/bin/
ln -s /usr/local/tools/fastqc/fastqc $HOME/bin/
ln -s /usr/local/tools/samtools/samtools $HOME/bin/
ln -s /Informatics_workshop/tools/annovar/table_annovar.pl $HOME/bin/
```

#2. Logout and relogin
1. `exit` # Exit the current ssh terminal
1. relogin

#3. tmux
"tmux" is a terminal multiplexer. It lets you switch easily between several programs in one terminal, detach them (they keep running in the background) and reattach them to a different terminal. And do a lot more.
1. `tmux` and enter
1. `binf`
1. `ctrl+b` and type `c` to open a new window 
1. `fastq`
1. `ctrl+b` and type `n` to swith to the next window
1. `ctrl+b` and type `p` to swith to the previous window
1. type `samt` and press `<tab>` key to see if `samtools` in `$PATH`

#4. File permission
1. `cd ~/`
1. `ls -la` or `l`
1. Examine the output more carefully. For example,
```bash
total 64K
drwx------. 14 hongc2 cc domain users 4.0K Jan 13 17:54 .
drwxr-xr-x. 37 root   root            4.0K Jan 11 17:53 ..
-rw-------.  1 hongc2 cc domain users  27K Jan 12 22:23 .bash_history
-rw-r-----.  1 hongc2 cc domain users   18 Sep  6 12:25 .bash_logout
-rw-r-----.  1 hongc2 cc domain users  247 Dec 29 12:27 .bash_profile
-rw-r-----.  1 hongc2 cc domain users  338 Dec 29 17:15 .bashrc
drwxrwx---.  2 hongc2 cc domain users  163 Jan 11 00:58 bin
drwxrwx---.  3 hongc2 cc domain users   18 Dec 26 21:23 .cache
drwxrwx---.  4 hongc2 cc domain users   30 Dec 26 23:18 .config
drwxrwx---.  3 hongc2 cc domain users   60 Dec 29 17:11 igv
drwxrwx---.  4 hongc2 cc domain users   37 Dec 29 17:10 .java
```
1. Q17: In the left column of the output, what is `-rwx`?
1. Group two people to work with
    1. Person A will share a file
    1. Person B will access/open the file
1. Person B: try to browse person A home directory
1. Person A: Share a file with your friend
    1. `cd /Informatics_workshop/shared_test`
    1. Create a file like `${USER}_shared.txt` (use `vi` editor)
1. Person B: try to access your neighbor shared file
    1. `cd /Informatics_workshop/shared_test`
    1. Can you access the file and open it?
    1. `less <neighbor>_shared.txt`
    
1. Person A: Change the file permission
    1. `chmod o-r ${USER}_shared.txt`
1. Person B:
    1. Switch to the next window in tmux (`ctrl+b`,`n`)
    1. Now, can you access the same file?
    1. `less <neighbor>_shared.txt`

#5. Identify your group ids
1. `groups`
1. `id -Gn`

#6. Bash script example
1. Q18: display first 10 read header for each fastq.gz file under the `fastq` directory
```bash
#!/bin/bash -l
for fastq_fn in $(ls ${HOME}/projects/bioinfo_2018/fastq/*.fastq.gz)
do
    zgrep -m10 '^@' ${fastq_fn}
done
```

#7. Grep example and pipe
1. Locate the VCF annotated table in example directory
    1. `binf`
    1. `l examples/annovar.chr22.snps.hg19_multianno.txt`
1. Build one line bash command satisfying all the following conditions
    1. variants where locate in exonic region
    1. nonsynonymous
    1. only columns 1, 2, 3, 4, 5, 6, 7, 10, and 17 are necessary
    1. gzip output file
```bash
grep -P '\texonic' annovar.chr22.snps.hg19_multianno.txt | grep -P '\tnonsynonymous' | cut -f 1,2,3,4,5,6,7,10,17 | gzip -fc > my_grep_test.gz 
```
#11. Access to Web resource
git, wget, and rsync
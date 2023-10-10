# customizations
# terminal command colored
export PS1="\[\033[36m\]\u\[\033[m\]@\[\033[32m\]\h:\[\033[33;1m\]\w\[\033[m\]\$ "
export CLICOLOR=1
export LSCOLORS=ExFxBxDxCxegedabagacad
alias ls='ls -G'
#alias diff='colordiff'
alias vi='vim'
alias grep='grep --color=auto'
alias fastjet-config='/Users/diptanilroy/ROOT_INSTALL/FASTJET/fastjet-3.3.4/fastjet-config'

# alias bnlRCF='ssh -L 10219:nx06.rcf.bnl.gov:22 -l droy1 rssh01.rhic.bnl.gov'
alias bnlRCF='ssh -L 10220:nx06.rcf.bnl.gov:22 -l droy1 ssh.sdcc.bnl.gov'
alias sPHENIXRCF='ssh -L 10220:nx06.rcf.bnl.gov:22 -l droy2 ssh.sdcc.bnl.gov'
alias hcal912='ssh droy2@hcal912gw.rhic.bnl.gov' #includes tunnel for VNC
# alias hcaldaq="ssh -M 20016 -L 5944:sphenixdaq:5901 sphenixdaq" ##includes tunnel for VNC
alias hcaldaq="ssh -L 5944:sphenixdaq:5901 sphenixdaq" ##includes tunnel for VNC

###########################################################################################
alias hometohome="ssh diptanilroy@192.168.1.124"
alias worldtohome="ssh diptanilroy@173.3.153.231"

####################### DAVIX #####################################
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/usr/local/lib/pkgconfig
# export OPENSSL_DIR = /usr/local/opt/openssl@1.1
# export OPENSSL_SUPPORT = -I$(OPENSSL_DIR)/include -L$(OPENSSL_DIR)/lib

######################## Change to working directory #################################
alias work="cd /Users/diptanilroy/Desktop/Y2020/STAR/STARAnalysis/2021"

alias D0Analysis="cd /Users/diptanilroy/Desktop/STAR-Workspace/D0Analysis"

######################## Load remote work directories ################################
alias rcfhome='umount -f /Users/diptanilroy/Desktop/RCF_MountPoint; sshfs droy1@sftp.sdcc.bnl.gov:/star/u/droy1/Y2019/STAR/ /Users/diptanilroy/Desktop/RCF_MountPoint'
alias gpfshome='umount -f /Users/diptanilroy/Desktop/GFPS_MountPoint; sshfs droy1@sftp.sdcc.bnl.gov:/gpfs01/star/pwg/droy1/ /Users/diptanilroy/Desktop/GFPS_MountPoint'

alias loadhome='umount -f /Users/diptanilroy/Desktop/HomePCAnalysisDirectory; sshfs diptanilroy@173.3.153.231:/Users/diptanilroy/Desktop/STAR-Workspace/ /Users/diptanilroy/Desktop/HomePCAnalysisDirectory'

alias sum="paste -sd+ - | bc"
# alias rcfhome= "umount -f /Users/diptanilroy/Desktop/RCF_MountPoint; sshfs droy1@sftp.sdcc.bnl.gov:/star/u/droy1/Y2019/STAR/ /Users/diptanilroy/Desktop/RCF_MountPoint"
# alias gpfshome= "umount -f /Users/diptanilroy/Desktop/GFPS_MountPoint; sshfs droy1@sftp.sdcc.bnl.gov:/gpfs01/star/pwg/droy1/ /Users/diptanilroy/Desktop/GFPS_MountPoint"

source /Users/diptanilroy/ROOT_INSTALL/ROOT/root-build/bin/thisroot.sh
# source /Applications/root_v6.18.04/bin/thisroot.sh
export PATH=/usr/local/bin:$PATH
export PATH=/Applications/CMake.app/Contents/bin:${PATH}
export PATH=/usr/bin:$PATH

export CC=/usr/bin/gcc
export CXX=/usr/bin/g++

# export FC=/Users/diptanilroy/ROOT_INSTALL/FORTRAN/usr/local

# export CC=/Applications/Xcode.app/Contents/Developer/usr/bin/gcc
# export CXX=/Applications/Xcode.app/Contents/Developer/usr/bin/g++

# export CC=/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/gcc
# export CXX=/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/gcc++
export PYTHIA8='/Users/diptanilroy/ROOT_INSTALL/PYTHIA/pythia8244'
export PYTHIA8DATA='/Users/diptanilroy/ROOT_INSTALL/PYTHIA/pythia8244/share/Pythia8/xmldoc'
export PATH=${PATH}:$PYTHIA8/bin

export FASTJET='/Users/diptanilroy/ROOT_INSTALL/FASTJET/fastjet-install'
#export FASTJET_CONTRIB='~/ROOT_INSTALL/FASTJET/fjcontrib-1.042'
export PATH=${PATH}:$FASTJET/bin

export LHAPDF6='/Users/diptanilroy/ROOT_INSTALL/LHAPDF/LHAPDF_install'
export LHAPDF6DATA='/Users/diptanilroy/ROOT_INSTALL/LHAPDF/LHAPDF_install/share/LHAPDF'
export PATH=${PATH}:$LHAPDF6/bin:$LHAPDF6DATA



export DYLD_LIBRARY_PATH=/usr/local/lib
export LD_LIBRARY_PATH=/usr/local/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$PYTHIA8/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PYTHIA8/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$FASTJET/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FASTJET/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$LHAPDF6/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LHAPDF6/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/diptanilroy/ROOT_INSTALL/ROOT/root-build/bin/thisroot.sh
#source /Users/diptanilroy/ROOT_INSTALL/ROOT/root-build/bin/thisroot.sh

save() {
	git add . && git commit -m "$1" && git push
}

pull(){
	git pull
}

### Ready-made function to copy files from BNL to local system
copyfromrcas () {

	sftp droy1@sftp.sdcc.bnl.gov <<< $'get -r '"$1"' '"$2"''

}

### Ready-made function to copy files to BNL from local system
copytorcas () {
	# scp -r droy1@rftpexp.rhic.bnl.gov:/star/u/droy1/Y2019/STAR/"$1" "$2"
	# scp -r $1 droy1@rftpexp.rhic.bnl.gov:$2
	#sftp droy1@sftp.rhic.bnl.gov <<< $'cd '"$2"'' <<< $'put -r '"$1"'' .

	sftp droy1@sftp.sdcc.bnl.gov:"$2" <<< $'put -r '"$1"''

}

homecopyfromhome(){
	scp -r diptanilroy@192.168.1.124:"$1" "$2"
}

homecopytohome(){
	scp -r $1 diptanilroy@192.168.1.124:$2
}

worldcopyfromhome(){
	scp -r diptanilroy@173.3.153.231:"$1" "$2"
}

worldcopytohome(){
	scp -r $1 diptanilroy@173.3.153.231:$2
}

### Ready-made function to copy files from BNL to local system
sphenixcopyfromrcas () {
    # scp -r droy1@rftpexp.rhic.bnl.gov:/star/u/droy1/Y2019/STAR/"$1" "$2"
    scp -r droy2@rftpexp.rhic.bnl.gov:"$1" "$2"

}

### Ready-made function to copy files to BNL from local system
sphenixcopytorcas () {
    # scp -r droy1@rftpexp.rhic.bnl.gov:/star/u/droy1/Y2019/STAR/"$1" "$2"
    scp -r "$1" droy2@rftpexp.rhic.bnl.gov:"$2"

}

# export DISPLAY=:10.0

# # >>> conda initialize >>>
# # !! Contents within this block are managed by 'conda init' !!
# __conda_setup="$('/Users/diptanilroy/opt/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
# if [ $? -eq 0 ]; then
#     eval "$__conda_setup"
# else
#     if [ -f "/Users/diptanilroy/opt/anaconda3/etc/profile.d/conda.sh" ]; then
#         . "/Users/diptanilroy/opt/anaconda3/etc/profile.d/conda.sh"
#     else
#         export PATH="/Users/diptanilroy/opt/anaconda3/bin:$PATH"
#     fi
# fi
# unset __conda_setup
# # <<< conda initialize <<<

# export PATH="/usr/local/bin:$PATH"
export PATH="/usr/local/Cellar/openssl@1.1/1.1.1i/bin:$PATH"
eval "$(thefuck --alias)"

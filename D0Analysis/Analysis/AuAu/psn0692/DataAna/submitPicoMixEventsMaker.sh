#########################################################################
# File Name: submitPicoD0AnaMaker.sh
# Created Time: Fri 08 May 2015 03:15:35 AM PDT
#########################################################################
#!/bin/bash

star-submit-template -template submitPicoMixEventsMaker.xml -entities listOfFiles=PicoDst.list

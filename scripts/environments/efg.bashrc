
# eFG Environment Config

#Edit this if required
SRC=~/src

#Launch efg in a bash subshell so we can exit without exiting login/launch shell

if [ $INIT_EFG ]; then
   . $SRC/ensembl-functgenomics/scripts/environments/efg.env

   #Reset INIT_EFG so any further subshell are not efg
   export INIT_EFG=
else
	#Lauch subshell with INIT_EFG set, and immediately reset in parent shell to
	#avoid further subshells launching efg
	alias efg='export INIT_EFG=1; bash -i; export INIT_EFG='
fi                      

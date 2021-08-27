#! /bin/bash
rm -f hostfile
touch hostfile
for i in csews{1..30}
do
	if ping -c1 -w1 $i &>/dev/null
	then 
		echo $i >> hostfile
	fi
done

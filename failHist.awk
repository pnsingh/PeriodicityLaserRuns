#!/usr/bin/awk -f
BEGIN{
	min=0
	max=710
	numSection=71
	count=0
	tot=0
}

{
	time[count,1]=$1
	time[count,2]=$2
	tot+=$2
	count++;

}

END{
	bin=int((max-min)/numSection)

	for(y=0;y<=10000;y++)
	{
		timeSection[y]=0

	}

	for(x=0;x<count;x++)
	{
			section=int((time[x,1]-min)/bin)	
			timeSection[section]+=time[x,2]
			print section*bin+min+bin "\t" timeSection[section]  > "hist.dat"

	}

}

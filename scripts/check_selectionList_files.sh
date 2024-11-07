# Takes the selectionList_*.txt as input

ip=$1
for i in `cat "$ip"`
do
	if [ ! -r $i ]; then
		echo "file $i doesn't exists"
	fi
done

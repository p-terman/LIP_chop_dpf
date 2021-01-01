#!/bin/bash

list=acq.list

dollar='$'
quote="'"
slash='\'
letter_r='r'

while read line
do
	name=${line%_*}
	echo $name	
	rm dummy_screen.sh
        echo -e "#!/bin/bash">> "dummy_screen.sh"
        echo -e "cd">> "dummy_screen.sh"
        echo -e "screen -d -m -S $name">> "dummy_screen.sh"
	echo -e "screen -S $name -X -p 0 stuff ${dollar}${quote}./Daemon_SRP.sh ${line} ${slash}${letter_r}${quote}">> "dummy_screen.sh"
	echo -e "screen -S $name -X -p 0 stuff $'exit ${slash}${letter_r}'">> "dummy_screen.sh"
	
	chmod +x dummy_screen.sh
	./dummy_screen.sh
	sleep 5

done < "${list}"



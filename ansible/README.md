To provision operating systems which only come with python3 installed (Ubuntu since 15.04),
execute the ansible command as following:

ansible-playbook site.yml -e 'ansible_python_interpreter=/usr/bin/python3'



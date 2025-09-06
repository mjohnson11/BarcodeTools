from sys import argv, stdout
import json

flist = argv[1:]

filenames = [f.split('/')[-1] for f in flist]

try:
    assert len(set(filenames)) == len(filenames)
except AssertionError:
    print('All filenames must be unique')


main_dict = dict()
for f in flist:
    with open(f, 'r') as infile:
        tmp = json.load(infile)
        main_dict[f.split('/')[-1]] = tmp
        
json_string = json.dumps(main_dict)
stdout.write(json_string)
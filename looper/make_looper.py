import glob
import numpy as np
import sys

img_string = sys.argv[1]
looper_filename = sys.argv[2]

print (('*' + model + '*' + station + '.png'))
files = np.sort(glob.glob(img_string))
print files
LIST = ""

for f in range(len(files)):
    string = 'window.temp_list[' + str(f) + ']= \"' + files[f] + '\";\n'
    LIST += string
print LIST
if looper_filename == 'looper_template.html':
    print "You don't want to overwrite the template, dude."
    sys.exit()
html_template = open('looper_template.html', 'r')
template = html_template.read()
template = template.replace("LIST", LIST)
template = template.replace("FIRSTIMAGE", files[0])
html_template.close()
out_html = open(looper_filename+'.html', 'w')
out_html.write(template)
out_html.close()


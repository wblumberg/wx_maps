import glob
import numpy as np
import sys
import os

img_string = sys.argv[1] + "*.png"
looper_filename = '/data/soundings/http/blumberg/' + sys.argv[2]
ROOT = 'http://sharp.weather.ou.edu/'

print "Searching for:", img_string

files = np.sort(glob.glob(img_string))
print "Found", len(files)

LIST = ""
for f in range(len(files)):
    link = files[f]
    abs_link = ROOT + link.split('http')[1]    
    string = 'window.temp_list[' + str(f) + ']= \"' + abs_link + '\";\n'
    LIST += string
print LIST

if looper_filename == 'looper_template.html':
    print "You don't want to overwrite the template, dude."
    sys.exit()

abs_path = os.path.dirname(os.path.abspath(__file__))
html_template = open(abs_path + '/looper_template.html', 'r')
template = html_template.read()
template = template.replace("LIST", LIST)
template = template.replace("FIRSTIMAGE", files[0])
template = template.replace("NIMGS", 'var imax = ' + str(len(files)) + ';')

html_template.close()

# Write the looper text out
print looper_filename + '.html'
out_html = open(looper_filename + '.html', 'w')
out_html.write(template)
out_html.close()


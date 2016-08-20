# Script that modifies python expressions in order to evaluate them faster
#
# USAGE: put the script into the home directory and run it. New folder
#        with modified expressions will be created
import re, os

directory_name = 'modified_expressions'

def read_in_pieces(file, piece_size = 1024*16):
    ''' reads file by chunks '''
    while True:
        piece = file.read(piece_size)
        if not piece:
            break
        yield piece

def modify_expression(fname):
    ''' modifies expression with given filename '''
    f = open('python_expressions/' + fname, 'r')
    t = open(directory_name + '/tmp', 'w')

    complex_values = {} # will store all complex literals in file
    gvars_to_add = '\n\n' # global (complex) variables declarations
    index_counter = 0 # index counter for complex variables

    # local variables declarations for 'exp()' expressions
    lvars_to_add = '\n\n    ex0 = exp(2 * mu[0])\n    ex1 = exp(2 * mu[1])\n    ex2 = exp(2 * mu[2])\n    ex3 = exp(2 * mu[3])\n'
    lvars_to_add += '    exm0 = exp(-2 * mu[0])\n    exm1 = exp(-2 * mu[1])\n    exm2 = exp(-2 * mu[2])\n    exm3 = exp(-2 * mu[3])\n'
    # 'power' variables declarations for powers of variables
    pvars_to_add = '\n\n'
    # variables that will have their powers simplified
    all_vars = ['r','K','k1','SP','SM','x[0]','x[1]','x[2]','zeta[0]','zeta[1]','zeta[2]','zeta[3]']
    for i in range(4):
        all_vars.append('ex' + str(i))
        all_vars.append('exm' + str(i))

    # stores the names of the 'power' variables
    pw_names = {}
    for var in all_vars:
        pw_names[var] = ['' for _ in range(6)]
        for p in range(2, 6):
            var_pw_name = 'p' + str(p) + re.sub('\]','',re.sub('\[','',var))
            pvars_to_add += '    ' + var_pw_name + ' = ' + var + '**' + str(p) + '\n'
            pw_names[var][p] = var_pw_name

    # read the file
    for piece in read_in_pieces(f):
        # find all complex literals
        for v in (re.findall(r"complex[(]\d+, \d+[)]", piece) + re.findall(r"complex[(]\d+, -\d+[)]", piece)):
            value = eval(v[7:])
            if not value in complex_values:
                new_var = 'c' + str(index_counter)
                index_counter += 1
                complex_values[value] = (new_var, v)
                gvars_to_add += new_var + ' = ' + v + '\n'

        # take the current piece and subtitute everything that needs to be substituted in it
        new_piece = piece

        # substitute complex literals
        for value, (var_name, v) in complex_values.items():
            new_piece = re.sub(re.sub('[(]', '[(]', re.sub('[)]', '[)]', v)), var_name, new_piece)

        # substitute 'exp()' expressions
        for i in range(4):
            new_piece = re.sub('exp[(]2 \* mu\[' + str(i) + '\][)]', 'ex' + str(i), new_piece)
            new_piece = re.sub('exp[(]-2 \* mu\[' + str(i) + '\][)]', 'exm' + str(i), new_piece)

        # substite powers
        for p in range(2, 6):
            for var in all_vars:
                new_piece = re.sub(re.sub('\[', '\\\[', re.sub('\]', '\\\]', var)) + ' \*\* ' + str(p), pw_names[var][p], new_piece)

        # write the substituted piece into temporary file
        t.write(new_piece)

    f.close()
    t.close()

    # now, it is time to add the variables declarations

    t = open(directory_name + '/tmp', 'r')
    targ = open(directory_name + '/' + fname, 'w') # target file

    # adds variables declarations in a hardcoded way
    first = True
    for piece in read_in_pieces(t):
        if first: # this is the first chunk
            lines = piece.split('\n')
            i = 0
            while lines[i][:4] == 'from': # add complex variables here
                i += 1
            j = i
            while lines[j][:3] != 'def': # add 'exp()' expressions here
                j += 1
            k = j
            while lines[k][:10] != '    DRMBAR': # add power variables here
                k += 1
            targ.write('\n'.join(lines[:i] + [gvars_to_add] + lines[i:j+1] + [lvars_to_add] + lines[j+1:k+1] + [pvars_to_add] + lines[k+1:]))
            first = False
        else:
            targ.write(piece)

    # close the files and clean up# this is the first chunk
    t.close()
    os.remove(directory_name + '/tmp')
    targ.close()

###################################################

# create directory if necessary
if not os.path.exists(directory_name):
    os.makedirs(directory_name)

# process all relevant expressions
expressions = [nm+str(i)+str(j)+str(k)+'.py' for nm in ['ddphis','dphis','ddgrams'] for i in range(1,4) for j in range(1,3) for k in range(1,3)]
expressions += ['dgrams'+str(i)+'.py' for i in range(1,4)] + ['NDD111.py']
for fname in expressions:
    print('Working on: ' + fname)
    modify_expression(fname)

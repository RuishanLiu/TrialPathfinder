import pandas as pd
import numpy as np
import datetime
import os

class cohort_selection(object):
    def __init__(self, patientids, name_PatientID='PatientID', lower=True):
        '''
        name_PatientID: the patientid name in all tables
        If lower=True, store all table names in lower case.
        '''
        patientids = list(patientids)
        patientids = pd.DataFrame(patientids, columns=[name_PatientID])
        patientids.index = patientids[name_PatientID]
        if patientids.shape[1] != 1 and len(patientids.shape)!=2:
            sys.exit('patientids must be one dimensional!')
            
        self.rules = dict()
        self.selects = dict()
        self.missing = dict()
        self.patientids = patientids
        self.lower = lower
        self.name_PatientID = name_PatientID
        self.tables = {}
        
    def self_name(self, name, lower=True):
        '''Given table name 'name', return the name stored in class self.name '''
        if lower:
            name = name.lower() # stored as lower class
        name = 'self.%s'%name
        return name

    def add_table(self, name, table):
        '''
        Add one table with name. 
        If lower=True, store all table names in lower case.
        Add a column of ordering = 0 - table.shape[0]-1
        '''
        if self.lower:
            name = name.lower()
        table = table.copy()
        table['ordering'] = np.arange(table.shape[0])
        table.index = np.arange(table.shape[0])
        self.tables[name] = table
        
    def add_tables(self, names, tables):
        '''
        Add tables with names (tolerant one table)
        names = [name1, name2, ...]
        tables = [table1, table2, ...]
        If lower=True, store all table names in lower case.
        '''
        if type(names) is str:
            self.add_table(names, tables, lower=self.lower)
        else:
            names = list(names)
            tables = list(names)
            if len(names) != len(tables):
                sys.exit('The numbers of names and tables do not match!')
            for name, table in zip(names, tables):
                self.add_table(name, table, lower=self.lower)
                
    def parse_line(self, line):
        '''
        Parse the line into logits=['|', '&'] and criteria=['table['re']==1', '...', '...']
        '''
        logits = []
        criteria = []

        # split the string to criteria
        starts = []
        ends = []
        i = 0
        starts.append(0)
        logits.append('')
        while i < len(line):
            if line[i:(i+2)] == 'OR':
                logits.append('|')
                ends.append(i)
                i += 2
                if i < len(line):
                    starts.append(i)
            elif line[i:(i+3)] == 'AND':
                logits.append('&')
                ends.append(i)
                i += 3
                if i < len(line):
                    starts.append(i)
            else:
                i += 1
        logits.append('')
        ends.append(len(line))

        # process the criteria
        while len(starts) > 0:
            start = starts.pop(0)
            end = ends.pop(0)
            crit = line[start:end]
            while (crit[0] == ' ') or (crit[-1] == ' ') or ((crit[0] == '(') and (crit[-1] == ')')):
                while (crit[0] == ' '):
                    crit = crit[1:]
                while (crit[-1] == ' '):
                    crit = crit[:-1]
                # remove bracket ( ) for each crit
                while (crit[0] == '(') and (crit[-1] == ')'):
                    crit = crit[1:]
                    crit = crit[:-1]
            criteria.append(crit)
                
        # Handle higher-level logics
        total_left = 0
        total_right = 0
        for i in range(len(criteria)):
            crit = criteria[i]
            n_left = np.sum([s=='(' for s in crit])
            n_right = np.sum([s==')' for s in crit]) 
            total_left += n_left
            total_right += n_right
            if n_left > n_right:
                for _ in range(n_left - n_right):
                    logits[i] = logits[i] + '('
                    criteria[i] = criteria[i][1:]
            elif n_left < n_right:
                for _ in range(n_right - n_left):
                    logits[i+1] = ')' + logits[i+1]
                    criteria[i] = criteria[i][:-1]
        return logits, criteria
        
    def retrieve_table(self, crit_string, return_index=False):
        '''
        Retrieve all table names in the crit_string as list. Also return crit_string with processed table name (self.lower)
        '''
        names_table = []
        starts = []
        ends = []
        i = len(crit_string)
        while i >=0:
            if crit_string[(i-2):i] == '[\'':
                end = i - 2
                j = end - 1
                while (crit_string[j].isalpha() or crit_string[j]=='_') and (j >= 0):
                    j = j - 1
                start = j + 1
                name_table = crit_string[start:end]
                if self.lower:
                    name_table = name_table.lower()
                names_table.append(name_table)
                i = start
                starts.append(start)
                ends.append(end)
            else:
                i = i - 1
        names_table = names_table[::-1] 
        starts = starts[::-1]
        ends = ends[::-1]
        if return_index:
            return names_table, starts, ends
        else:
            return names_table
        
    def add_prefix(self, name, table):
        '''
        Add prefix name_ to every column name of table, except for self.name_PatientID
        '''
        columns = table.columns
        columns_processed = []
        for col in columns:
            if col != self.name_PatientID:
                col = '%s_%s' % (name, col)
            columns_processed.append(col)
        table.columns = columns_processed
        return table
    
    def crit_to_command(self, crit, name_exec='self.exec_tables', single=True):
        '''
        Turn input crit to command could be exec
        If single = True, only one table is in the criterion
        '''
        def extract_date(j, non_table):
            '''
            extract the date outof non_table starting from j.
            Return new_j (after ')'), string of dates (e.g. 50)
            j: non_table[j] == (
            '''
            if non_table[j] != '(':
                sys.exit('Format for date must be @YEARS(4) or @MONTHS(5) or @DAYS(30)!')
            j += 1
            string_date = ''
            while non_table[j] != ')' and j < len(non_table):
                string_date += non_table[j]
                j += 1
            if non_table[j] != ')':
                sys.exit('Format for date must be @YEARS(4) or @MONTHS(5) or @DAYS(30)!')
            j += 1 # new start
            try:
                string_date = int(string_date)
            except:
                sys.exit('Format for date must be @YEARS(4) or @MONTHS(5) or @DAYS(30)!')
            return j, string_date
                
            
        names_table, starts, ends = self.retrieve_table(crit, return_index=True)
        command = crit[:starts[0]]
        sign_date = ('@YEARS' in crit) or ('@MONTHS' in crit) or ('@DAYS' in crit)
        for i in range(len(starts)):
            name_table = names_table[i]
            start_next = (len(crit)+1) if i==(len(starts)-1) else starts[i+1]
            # table part
            if single:
                name_table = '%s[\'%s\']' % (name_exec, name_table)
                end_pre = ends[i]
            else:
                name_table = '%s[\'table_merge\'][\'%s_' % (name_exec, name_table)
                end_pre = ends[i]+2
            non_table = crit[end_pre:start_next]
            
            # includes date. @YEARS()  @MONTHS()  @DAYS()
            if sign_date:
                # first find the table['feature']
                j = 0
                while non_table[j:(j+2)] != '\']' and j < len(non_table):
                    j += 1
                if j == len(non_table) - 1:
                    sys.exit('Format is wrong - table[\'feature\'].')
                name_table += non_table[:(j+2)]
                non_table = non_table[(j+2):]
                name_table = 'pd.to_datetime(%s)' % name_table
                non_table_date = ''
                j = 0
                while j < len(non_table):
                    multiply = 0
                    if non_table[j:(j+6)] == '@YEARS':
                        j += 6
                        multiply = 365
                    elif non_table[j:(j+7)] == '@MONTHS':
                        j += 7
                        multiply = 30
                    elif non_table[j:(j+5)] == '@DAYS':
                        j += 5
                        multiply = 1
                    if multiply > 0:
                        j, date = extract_date(j, non_table)
                        non_table_date += 'datetime.timedelta(%d)' % (multiply * date)
                    else:
                        non_table_date += non_table[j]
                        j += 1 
                non_table = non_table_date
      
            command +=  name_table
            command += non_table 
            
        #### if MIN or MAX #####
        comp = None
        if command[:3] == 'MIN':
            comp = 'min'
        elif  command[:3] == 'MAX':
            comp = 'max'
        if comp is not None:
            command = command[3:]
            
        # ABS()
        command = command.replace('ABS', 'abs')
        
        return command, comp
    
    
    def process_minmax(self, comp, select):
        '''
        Process when MIN() MAX() exists.
        select: values to be operated by argmin or argmax
        Patientwise. Return whether selected in self.exec_tables['table_merge']'s patientid
        '''
        values = select.copy()
        select.iloc[:] = False
        patientid_set = list(set(list(self.exec_tables['table_merge'].loc[values.index, self.name_PatientID])))
        for patientid in patientid_set:
            table_idx = list(self.exec_tables['table_merge'].loc[self.exec_tables['table_merge'][self.name_PatientID]==patientid].index)
            if len(table_idx) > 0:
                values_sub = list(values.loc[table_idx])
                if comp == 'min':
                    ind = np.argmin(values_sub)
                else:
                    ind = np.argmax(values_sub)
                idx = table_idx[ind]
                select[idx] = True
        return select 

    def execute_rule(self, rules, force=False):
        '''
        Execute rule and return whether the patient in patientsid selected or not
        If force=True, overwrite the rule in self.rules if the same name exists.
        '''
        self.exec_tables = {}
        name_exec = 'self.exec_tables'
        self.exec_var = {}
        
        rules = rules.split('\n#')
        if len(rules) == 1:
            print('Error: rule must start with #. Rule not added!')
            return None
        
        select_result = None
        name = None
        for rule in rules:
            if len(rule) == 0:
                continue
            rule = rule.split('\n')
            if name is None:
                rule = [r for r in rule if r!='']
                name = rule[0]
                if name in self.rules:
                    if force:
                        print('Warning: Rule named "%s" already exists. Now overwrite (force=True)' % name)
                    else:
                        print('Warning: Rule named "%s" already exist. Rule not added! Set force=True or remove this rule first!' % name)
                        return name, self.selects[name], self.missing[name]
                continue
            ie = rule[0].lower()
            flag_in_missing = True
            if 'missing exclude' in ie:
                flag_in_missing = False
            ie = ie.replace(' ', '')
            ie = ie[:2]
            if ie not in ['in', 'ex']:
                print('Error: #Inclusion, #Exclusion only. Rule not added!')
                return None 
            if len(rule) < 2:
                print('Error: rule format is wrong. Rule not added!')
                return None 
        
            #### Execute criteria ####
            lines_raw = rule[1:]
            lines = []
            for line in lines_raw:
                if line != '':
                    line = line.replace('’', '\'')
                    line = line.replace('‘', '\'')
                    lines.append(line)

            # Initiate All tables in the criteria
            crits = ' '.join(lines)
            names_table = self.retrieve_table(crits)
            names_table_set = list(set(names_table))
            for name_table in names_table_set:
                self.exec_tables[name_table] = self.tables[name_table]

            last = False # whether the last line
            single = True # whether only one is per criterion
            merge = False # wheter merged table or not
            selects_final = []
            for i, line in enumerate(lines):
                if i == len(lines) - 1:
                    last = True
                logits, criteria = self.parse_line(line)
                selects = []
                for crit in criteria:
                    names_table = self.retrieve_table(crit)
                    name_table = names_table[0] 
                    # Check whether to start single mode
                    if len(names_table) > 1:
                        single = False
                    if (not single) and (not merge):
                        tables_use = [self.add_prefix(name, self.exec_tables[name].copy()) for name in names_table_set]
                        table_merge = tables_use[0]
                        for table in tables_use[1:]:
                            table_merge = pd.merge(table_merge, table, on=self.name_PatientID).copy()
                        self.exec_tables['table_merge'] = table_merge.copy()
                        merge = True

                    # Only one table appears = select: whether selected for name_table (the first from left to right) [True, False, ...]
                    # More tables appears = select: whether selected for the merged table [True, False, ...]
                    
                    command, comp = self.crit_to_command(crit, name_exec=name_exec, single=single)  
                    exec_string = 'self.exec_var[\'select\'] = %s' % (command)
                    exec(exec_string)
                    select = self.exec_var['select']
                    if comp is not None:
                        select = self.process_minmax(comp, select)
                    
                    if not single:
                        self.exec_tables['table_merge'] = self.exec_tables['table_merge'].loc[select]
                        orders_in = self.exec_tables['table_merge'].loc[select, '%s_ordering' % name_table]
                        # select: whether selected for name_table
                        select = self.exec_tables[name_table]['ordering'].isin(orders_in)

                    if last:
                        # select = whether selected for patientids
                        exec_string = 'self.exec_var[\'select\'] = self.exec_tables[\'%s\'].loc[select, \'%s\']' % (name_table, self.name_PatientID)
                        exec(exec_string)
                        select = self.exec_var['select']
                        select = self.patientids[self.name_PatientID].isin(select)
                        selects.append(select)
                    else:
                        # select = whether selected for name_table
                        selects.append(select)

                # Process with Logits. select_final = select1 & select2 | ...
                exec_string_logit = 'self.exec_var[\'select_final\'] = '
                for i_select, select in enumerate(selects):
                    exec_string_logit += logits[i_select]
                    name_select = 'self.exec_var[\'select_%d\']' % (i_select)
                    exec_string = '%s = select' % (name_select)
                    exec(exec_string)
                    exec_string_logit += name_select
                exec_string_logit += logits[-1]
                exec(exec_string_logit)
                select_final = self.exec_var['select_final']
                
                if last:
                    exec_string = 'self.exec_var[\'missing\'] = ~self.patientids[self.name_PatientID].isin(self.exec_tables[\'%s\'][self.name_PatientID])' % (name_table)
                    exec(exec_string)
                    missing = self.exec_var['missing']
                    # Patient missing in any table of the last statement
                    for name_table_last in names_table[1:]:
                        exec_string = 'self.exec_var[\'missing\'] = ~self.patientids[self.name_PatientID].isin(self.exec_tables[\'%s\'][self.name_PatientID])' % (name_table_last)
                        exec(exec_string)
                        missing = missing | self.exec_var['missing']

                # Operation on name_table (the first from left to right)
                if not last:
                    # operation on name_table. subset by select_final
                    exec_string = 'self.exec_tables[\'%s\'] = self.exec_tables[\'%s\'].loc[select_final]' % (name_table, name_table)
                    exec(exec_string)
                    # get the patientid from new name_table
                    exec_string = 'self.exec_var[\'select_final\'] = self.exec_tables[\'%s\'][\'%s\']' % (name_table, self.name_PatientID)
                    exec(exec_string)
                    # whether self.patientid in patientid
                    select_final = self.exec_var['select_final']
                    select_final = self.patientids[self.name_PatientID].isin(select_final)

                selects_final.append(select_final)

            select_all = selects_final[0]
            for select in selects_final[1:]:
                select_all = select_all & select

            if ie == 'ex':
                select_all = ~select_all
        
            # Merge the results for all 
            if select_result is None:
                select_result = select_all
            else:
                select_result = select_result & select_all
                
            if flag_in_missing:
                select_result |= missing
        
        if select_result is None: # select not generated
            return None
        
        return name, select_result, missing
           
    def add_rule(self, rule, force=False):
        '''
        Add rule. self.rules[rule] = select (patientids whether selected)
        If force=True, overwrite the rule in self.rules if the same name exists.
        '''
        result = self.execute_rule(rule, force=force)
        if result is not None:
            name, select, missing = result
            self.rules[name] = rule
            self.selects[name] = select
            self.missing[name] = missing
            return name, select, missing
        else:
            return None, None, None
        
        
    def add_rules(self, rules, force=False):
        '''
        Add rules
        '''
        if type(rules) is str:
            self.add_rule(rules, force=force)
        else:
            for rule in rules:
                self.add_rule(rule, force=force)
            
    def selection(self, name_rules):
        '''
        Given list of name_rules, return whether self.patients are selected or not.
        '''
        if type(name_rules) is str:
            select = self.selects[name_rules].copy()
        else:
            select = self.selects[name_rules[0]].copy()
            for name_rule in name_rules[1:]:
                select &= self.selects[name_rule]
        select_id = list(self.patientids.loc[select, self.name_PatientID])
        self.select_id = select_id
        return select_id
    
    
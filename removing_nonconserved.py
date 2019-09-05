#!/usr/bin/env python
"ordering as well as rotation of the genomes is done for almost conserved genes"
"missing genes are stored in sorted order"
from xlrd import open_workbook
import xlsxwriter
    
    
wb = open_workbook("FILE.xlsx")
workbookfinal = xlsxwriter.Workbook("removed_not_conserved"+'.xlsx')
worksheetfinal = workbookfinal.add_worksheet()


values=[]
counter=0

for s in wb.sheets():
    for row in range(s.nrows):
        counter=0
        row_value = []
        for col in range(s.ncols):
            value  = (s.cell(row,col).value)
            if value=='-':
                counter=counter+1
            try : value = int(value)
            except : pass
            row_value.append(value)

        if counter==1 or counter==0:
            values.append(row_value)

            
        
for i in range(0,len(values)):
    val=values[i]
    for j in range(0,len(val)):
        worksheetfinal.write(i,j,val[j])

workbookfinal.close()




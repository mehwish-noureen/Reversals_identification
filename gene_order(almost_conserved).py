#!/usr/bin/env python
"ordering as well as rotation of the genomes is done for almost conserved genes"
"missing genes are stored in sorted order"
from xlrd import open_workbook
import xlsxwriter



def rotate(genome,value):
    ind1=genome.index(1)
    ind2=genome.index(value)
    if ((ind1!=len(genome)-1 and ind2!=0) and (ind1<ind2)):
        part1=genome[0:ind1+1]
        part2=genome[ind1+1:len(genome)]
        part1.reverse()
        part2.reverse()
        genome=part1+part2
    elif ((ind1!=len(genome)-1 and ind2!=0) and (ind1>ind2)):
        part1=genome[ind1:len(genome)]
        part2=genome[0:ind2+1]
        genome=part1+part2
    elif ind1==0 and ind2!=len(genome)-1:
        part1=genome[0:ind2-1]
        part2=genome[ind2:len(genome)-1]
        part2.reverse()
        genome=part1+part2
    else:
        genome.reverse()

    return genome
    
    
wb = open_workbook("Genomic_Position.xlsx")
workbookfinal = xlsxwriter.Workbook("order_new_rotated"+'.xlsx')
worksheetfinal = workbookfinal.add_worksheet()
worksheetfinal1 = workbookfinal.add_worksheet()

values=[]
counter=[]

for s in wb.sheets():
    for col in range(s.ncols):
        col_value = []
        for row in range(s.nrows):
            value  = (s.cell(row,col).value)
            try : value = int(value)
            except : pass
            col_value.append(value)

        values.append(col_value)
        if col!=0:
            dd=col_value.count('-')
            counter.append(dd)


Gene_order=[]
order=[]


for z in range(len(values)-1):
    c1=[values[0]]+[values[z+1]]
    new_sort=[valx for (valy,valx) in sorted(zip(c1[1],c1[0]))]
    
    if z==0:
        order=range(1,len(new_sort)+1)
        order=[order]+[new_sort]
        g=order[0]
        Gene_order=[g]

    else:
        for r in range(0,len(new_sort)):
            num=new_sort[r]
            ind=order[1].index(num)
            ll=order[0]
            num1=ll[ind]
            new_sort[r]=num1

        Gene_order.append(new_sort)


xx=len(Gene_order)
last_val=len(values[1])-counter[0]
new_order=[]
missing_val=[]

for g in range(len(Gene_order)):
    kk=Gene_order[g]
    kk_new=kk[0:len(kk)-counter[g]]
    miss=kk[(len(kk)-counter[g]):len(kk)]
    new_order.append(kk_new)
    if miss==[]:
        miss=[0]
    missing_val.append(sorted(miss))

for ii in range(0,xx):
    val_order=new_order[ii]
    if val_order[0]!=1 or val_order[len(val_order)-1]!=last_val:
        print "hello"
        r_gene=rotate(val_order,last_val)
        new_order[ii]=r_gene
        worksheetfinal.write_column(0,ii,r_gene)
        worksheetfinal1.write_column(0,ii,missing_val[ii])
    else:    
        worksheetfinal.write_column(0,ii,val_order)
        worksheetfinal1.write_column(0,ii,missing_val[ii])

workbookfinal.close()




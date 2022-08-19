function showDiv(chr){
	line=fileText.split(/\n/);

	rows=line.length;
	cols=16;
	// alert(rows+'\n'+cols)
	var tab='<table border=1 width=700'
	tab+='<tr>';
	tab+="<td style='background:#6b778d'>chromosome</td>";
	tab+="<td style='background:#6b778d'>Pos</td>";
    tab+="<td style='background:#6b778d'>rsID</td>";
	tab+="<td style='background:#6b778d'>Ref</td>";
	tab+="<td style='background:#6b778d'>Alt</td>";
	tab+="<td style='background:#6b778d'>Ref Count</td>";
	tab+="<td style='background:#6b778d'>Alt Count</td>";
	tab+="<td style='background:#6b778d'>p.Val</td>";
    tab+="<td style='background:#6b778d'>Region</td>";
	tab+="<td style='background:#6b778d'>cCRE</td>";
	tab+="<td style='background:#6b778d'>Motif Overlap</td>";
	tab+="<td style='background:#6b778d'>GenomeBroswer link</td>";
	tab+="<td style='background:#6b778d'>VariantViewer link</td>";
	tab+="<td style='background:#6b778d'>GTEx</td>";
	tab+="<td style='background:#6b778d'>GWAS</td>";
	tab+='</tr>';

	for(var i=0;i<rows;i++){
		line1=line[i];
		chrom=line1.split(/ /)[0];
		if (chrom==("chr"+chr.innerHTML)){
			pos=line1.split(/ /)[1];
            rsID=line1.split(/ /)[2];
			ref=line1.split(/ /)[3];
			alt=line1.split(/ /)[4];
			ref_count=line1.split(/ /)[5];
			alt_count=line1.split(/ /)[6];
			pvalue=line1.split(/ /)[11];
            ann=line1.split(/ /)[14];
			ccre=line1.split(/ /)[15];
			ccre_id=line1.split(/ /)[16];
			motif=line1.split(/ /)[17];
			gb=line1.split(/ /)[20];
			vv=line1.split(/ /)[21];
			gtex=line1.split(/ /)[18];
			gwas=line1.split(/ /)[19];

			p_str=pvalue.split(/e/);
			p=parseFloat(p_str[0]);
			p=p.toFixed(4)
			if (p_str.length==2){
				p=p+'e'+p_str[1]
			}

			arr=[chrom, pos, rsID, ref,alt, ref_count, alt_count, p, ann, ccre, ccre_id, motif, 'Genome Broswer', 'Variant Viewer', gtex, gwas]
			tab+='<tr>'
		    for(var j=0;j<cols;j++){
		    	if (j!=10){
			    	if (j<12||j>13){
			    		if (j==2 & arr[j]!='-'){
			    			tab+="<td style='background:#6b778d'><a href=https://www.ncbi.nlm.nih.gov/snp/"+arr[j]+">"+arr[j]+"</a></td>"
			    		}
			    		else if (j==9 & arr[j]!='-'){
			    			tab+="<td style='background:#6b778d'><a href=https://screen.encodeproject.org/search/?q="+arr[j+1]+"&assembly=GRCh38&uuid=5d366e00-5f48-4da2-86cb-77c333fb5651>"+arr[j]+"</a></td>"
			    		}
			    		else{
			    			tab+="<td style='background:#6b778d'>"+arr[j]+"</td>"
			    		}
			    	}
			    	else if (j==12){
			    		tab+="<td style='background:#6b778d'><a href="+gb+">"+arr[j]+"</a></td>"
			    	}
			    	else{
			    		tab+="<td style='background:#6b778d'><a href="+vv+">"+arr[j]+"</a></td>"
			    	}
				}
			}
			tab+='</tr>'
			
		}
	}
	tab+='</table>';

	var showBox=document.getElementById("showBox");
	showBox.innerHTML=tab;
}


#include<stdio.h>
#include<math.h>
#include<stdlib.h>
main()
{
int i,cnt,file_mode,no_line,c;
float mjd,obsL,freq,ch_wid;
double tsamp;
int obsL1,bits;
int mjd1,bitshift;
char tmp[100],ch[50];
char date[50],utc[50],ra[50],dec[20],src[50],time_res[20],freq_hdr[20],junk[50];
char* p;
FILE *fp,*fp1,*out;


  fp=fopen("template.gmrt_hdr","r");
  fp1=fopen("param","r");
  out=fopen("gmrt_hdr","w");

fscanf(fp1,"%s%d%s%s%s%s%s%s%s",date,&mjd1,junk,utc,src,ra,dec,time_res,freq_hdr);
printf("%s %d %s %s %s %s %s %s\n",date,mjd1,utc,src,ra,dec,time_res,freq_hdr);

//mjd1 = (int)(mjd);

  if (fp== NULL) {
	          printf("File 'template.gmrt_hdr' failed to open\n");
		  exit (1);
                  } 
  if (out== NULL) {
	          printf("File 'gmrt_hdr' failed to write\n");
		  exit (1);
                  } 
i=0;
do{
  p=fgets(tmp,100,fp);
	if (p != NULL) {
cnt=0;
	switch(i)
	{
		case 6 :
			 do{
			    fprintf(out,"%c",tmp[cnt]); 
			    cnt++;
			    }while(tmp[cnt]!='\n');
                   fprintf(out,"%s \n",date);
		   break;
		case 11 :
                         do{
                            fprintf(out,"%c",tmp[cnt]);
                            cnt++;
                            }while(tmp[cnt]!='\n');
                   fprintf(out,"%s \n",freq_hdr);
		   break;
		case 12 :
                         do{
                            fprintf(out,"%c",tmp[cnt]);
                            cnt++;
                            }while(tmp[cnt]!='\n');
                   fprintf(out,"%s \n",time_res);
                   break;   
		case 16 :
			 do{
			    fprintf(out,"%c",tmp[cnt]); 
			    cnt++;
			    }while(tmp[cnt]!='\n');
                   fprintf(out," %d \n",mjd1);
                   fprintf(stdout," %d %f \n",(int)mjd,mjd);
		   break;
		case 17 :
			 do{
			    fprintf(out,"%c",tmp[cnt]); 
			    cnt++;
			    }while(tmp[cnt]!='\n');
                   fprintf(out," %s \n",utc);
		   break;
		case 18 :
			 do{
			    fprintf(out,"%c",tmp[cnt]); 
			    cnt++;
			    }while(tmp[cnt]!='\n');
                   fprintf(out,"%s \n",src);
		   break;
		case 19 :
			 do{
			    fprintf(out,"%c",tmp[cnt]); 
			    cnt++;
			    }while(tmp[cnt]!='\n');
                   fprintf(out," %s, %s \n",ra,dec);
		   break;
		default :
                   fprintf(out,"%s",tmp);
                   break;
	}
    }
   i++;
   } while (p != NULL);

    printf("Done Creating gmrt_hdr file!!\n");
    fclose(fp);
    fclose(out);

}

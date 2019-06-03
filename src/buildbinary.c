/*==================================================================================================
	C component of build.binary.R 

	This C function fills a matrix considering dependancy between sites.
	It construct the binary matrix containing the presence(1) or the absence (0) of a link to a site
====================================================================================================*/

/*================
	C library used
==================*/
# include <stddef.h>
# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <assert.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
# include <R_ext/Utils.h>
# include "adesub.h"


/*========================================
	Term declaration - object created in R
==========================================*/

void buildbinary ( int *nrowlinkR,  int *linkR,
				   int *pointsorderR,  int *lengthpoR,
				   int *nsiteR,  int *matR);




/*======================================
	Term definition - object created in R
========================================*/
void buildbinary ( int *nrowlinkR,  int *linkR,
				   int *pointsorderR,  int *lengthpoR,
				   int *nsiteR,  int *matR)

{

/*=======================================
	Term definition - object created in C
=========================================*/

 int a,b,c,d,e,f,g,h,i,sizetokeep,*tmp,sizetmp,sizetmpold,*tmp2,sizetmp2,**linktmptmp,*line,*linetmp,
			 **linktmp,nsiteC,**linkC,nrowlinkC,*pointsorderC,lengthpoC,**matC;


	/*=====================
		Pass R objects to C
	=======================*/

	/*-------------------------------
		Assign nrowlinkR to nrowlinkC
	---------------------------------*/
	nrowlinkC = *nrowlinkR;

	/*-------------------------
		Build linkC
		Assign linkC to linktmp
	---------------------------*/

	tabintalloc(&linkC, nrowlinkC, 4);
	tabintalloc(&linktmp, nrowlinkC, 4);

	a = 0;
	for (b = 1; b <= 4; b++) {
		for (c = 1; c <= nrowlinkC; c++) {
			linkC[c][b] = linkR[a];
			linktmp[c][b] = linkR[a];
			a = a + 1;
		}
	}
		
	/*-------------------------------
		Assign lengthpoR to lengthpoC
	---------------------------------*/
	lengthpoC = *lengthpoR;
	
	/*----------------------
		Rebuild pointsorderR
	------------------------*/
	vecintalloc(&pointsorderC,lengthpoC);
	
	a = 0;
	for (b = 1; b <= lengthpoC; b++) {
		pointsorderC[b] = pointsorderR[a];
		a = a + 1;
	}
	
	/*-------------------------
		Assign nsiteR to nsiteC
	---------------------------*/
	nsiteC = *nsiteR;
	
	/*------------
		build matC
	--------------*/
	tabintalloc(&matC, nsiteC,nrowlinkC);

	/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
		Building the binary matrix (actual start of the real procedure)
	-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
	
	vecintalloc(&tmp,1);

	/*Loop selecting each link*/
	for (a = 1; a <= nrowlinkC; a++){
		
		/*## Fill up matC considering each link and the dependency of each links ##*/
		/*## (If one site is related to another via a middle one (1-> 2-> 3) is considered a dependancy) ##*/
		/*## This type of loop for situations where sites aren't necessarily ordered nicely ##*/
		/*Loop starting for each sites*/
		b = 1;
		while(b <= lengthpoC){
			c = b;
			b = pointsorderC[c];


			sizetokeep = nrowlinkC;
			
			/*------------
				build line
			--------------*/
			vecintalloc(&line,nrowlinkC);
			
			for (d = 1; d <= nrowlinkC; d++){
				line[d] = linktmp[d][1];
			}
			

			/*Initiate sizetmp and sizetmpold*/
			sizetmp = 1;
			sizetmpold = 1;
				
			/*Select the links which have the site number "b" in the loop*/
			tmp[1] = a;
			
			/*## Loop that goes on until sizetmp equal 0 ##*/
			/*## (Until linktmp is empty) ##*/
			while(sizetmp != 0){
				d = 1;
				while(d <= sizetmp){
					e = d;
					d = tmp[e];
					/*Select the object at the end of the link considered*/
					f = linktmp[d][3];
					
					/*Select the link considered*/
					g = linktmp[d][1];
					
					/*Put 1 in the cell of matC where the object at the end of the link and the link meets*/
					matC[f][g] = 1;
					
					/* Find which sites at the end of the link considered is the same as the starting link*/
					vecintalloc(&tmp2,nrowlinkC);
					sizetmp2 = 0;
					h = 1;
					for(i = 1; i <= sizetokeep; i++){
						if(linktmp[i][4]==f){
							tmp2[h] = line[i];
							h = h + 1;
							sizetmp2 = sizetmp2 + 1;
						}
					}
					
					/*## if tmp2 is of length not nul, modify linktmp ##*/
					/*## fourth column with b and first column with g ##*/
					if(sizetmp2 > 0){
						for(h = 1; h <= sizetmp2; h++) {
							linktmp[tmp2[h]][4] = b;
							linktmp[tmp2[h]][1] = g;
						}
					}
					d = e;
					d = d + 1;
					freeintvec(tmp2);
				}
				
				/* Copy linktmp in linktmptmp while removing object used in tmp */
				tabintalloc(&linktmptmp, sizetokeep-sizetmp, 4);
				
				for(e = 1; e <= 4; e++){
					g = 1;
					for(d = 1; d <= sizetokeep; d++){
						/* This little twist here is for situations where tmp is larger than 1 */
						i = 0;
						for(h = 1; h <= sizetmp; h++){
							if(d != tmp[h]){
								i = i + 1;
							}
						}
						if(i == sizetmp){
							linktmptmp[g][e] = linktmp[d][e];
							g = g + 1;
						}
					}
				}

				/* Reconstruct line with a proper size */
				vecintalloc(&linetmp,sizetokeep-sizetmp);

				for(d = 1; d <= sizetokeep-sizetmp; d++){
					linetmp[d] = line[d];
				}
				
				freeintvec(line);
				vecintalloc(&line,sizetokeep-sizetmp);
				
				for(d = 1; d <= sizetokeep-sizetmp; d++){
					line[d] = linetmp[d];
				}
				
				freeintvec(linetmp);
				
				/* Reconstruct linktmp with a proper size */
				freeinttab(linktmp);
				tabintalloc(&linktmp, sizetokeep-sizetmp, 4);
				
				for(d = 1; d <= sizetokeep-sizetmp; d++){
					for(e = 1; e <= 4; e++){
						linktmp[d][e] = linktmptmp[d][e];
					}
				}
				
				freeinttab(linktmptmp);
				
				/* Reconstruct tmp */
				freeintvec(tmp);
				vecintalloc(&tmp,sizetokeep-sizetmp);

				e = 1;
				for(d = 1; d <= sizetokeep-sizetmp; d++){
					if(b == linktmp[d][4]){
						tmp[e] = line[d];
						e = e + 1;
					}
				}
				
				/* Reconstruct sizetmp, sizetmpold and sizetokeep */
				sizetmpold = sizetmp;
				sizetmp = e - 1;

				sizetokeep = sizetokeep - sizetmpold;
				
				/* If sizetokeep is the of the same length as sizetmp, stop the loop */
				if(sizetokeep == sizetmp){
					break;
				}
			}

			/* Reconstruct linktmp from scratch... ok well... from linkR  */
			freeinttab(linktmp);
			tabintalloc(&linktmp, nrowlinkC, 4);

			d = 0;
			for(e = 1; e <= 4; e++) {
				for (f = 1; f <= nrowlinkC; f++) {
					linktmp[f][e] = linkR[d];
					d = d + 1;
				}
			}
			
			freeintvec(line);
			
			b = c;
			b = b + 1;
			
		}
	}

	/*------------
		build matR
	--------------*/
	
	a = 0;
	for (b = 1; b <= nsiteC; b++) {
		for (c = 1; c <= nrowlinkC; c++) {
			matR[a] = matC[b][c];
			a = a + 1;
		}
	}

	/* Free memory taken by object matC */
	freeinttab(matC);
	freeinttab(linkC);
	freeinttab(linktmp);
	freeintvec(tmp);
}

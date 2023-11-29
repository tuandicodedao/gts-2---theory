#include<stdio.h>
#include<conio.h>
#include<math.h>
#define M 100
#define E 0.0001
float kiemtra(float b[M][M],int h1,int h2, int k)
{
int x;float temp;
	temp=fabs(b[1][h1]-b[1][h2]);
	for (x=1;x<=k;x++)
		if (temp<fabs(b[x][h1]-b[x][h2])) temp=fabs(b[x][h1]-b[x][h2]);
	return temp;
}
int main()
{
//khai bao bien
int i,j,k,n,th,tam;
float a[M][M], b[M][M],s,x,y,z,temp;
FILE *f1, *f2;
//nhap du lieu
f1=fopen("input.txt","r");
fscanf(f1,"%d",&n);
printf("%d\n",n);
for (i=1;i<=n;i++)	
	for (j=1;j<=n;j++) fscanf(f1,"%f ",&a[i][j]);
fclose(f1);	
for (i=1;i<=n;i++){
	for (j=1;j<=n;j++) printf("%3.2f  ",a[i][j]);
	printf("\n");}
printf("\n");
f2=fopen("output.txt","w");
for (int lap=1;lap<=n;lap++){
	for (i=1;i<=n;i++) b[i][0]=1;        //khoi tao vector Y
	for (k=1;k<=3;k++)			//tinh 3 vector dau trong day
		for (i=1;i<=n;i++){
			s=0;
			for (j=1;j<=n;j++)
				s+=a[i][j]*b[j][k-1];
			b[i][k]=s;}
	j=3;
	for(;;){
		j++;
		for (i=1;i<=n;i++){	//tinh cac vector tiep theo
			s=0;
			for (k=1;k<=n;k++)
				s+=a[i][k]*b[k][j-1]; b[i][j]=s;}
		temp=b[1][j];
		for (i=1;i<=n;i++)if(fabs(b[i][j])>fabs(temp)) temp=b[i][j];
		for (i=1;i<=n;i++) b[i][j]=b[i][j]/temp;
		if (kiemtra(b,j,j-1,n)<=E){th=1;break;}
		if (kiemtra(b,j,j-2,n)<=E&&(b,j,j-1,n)>1){th=3;break;}
		if (kiemtra(b,j,j-1,n)-kiemtra(b,j-1,j-2,n)>0&&(kiemtra(b,j,j-2,n)-kiemtra(b,j-1,j-3,n)>0)&&(j>7))
			{th=4;break;}
		if (j>=100){th=4;break;}}
	printf("\nday vector tinh duoc\n");
	for (i=1;i<=n;i++){
		for (k=1;k<=j;k++) printf("%10f ",b[i][k]);
		printf("\n");}
	switch(th){
		case 1:			
		for (i=1;i<=n;i++)
			{s=0;
			for (k=1;k<=n;k++)
				s+=a[i][k]*b[k][j-1];
				b[i][j]=s;}
		temp=b[1][j];tam=1;
		for (i=1;i<=n;i++) if (b[i][j]>temp) {temp=b[i][j];tam=i;}
		printf("gia tri rieng thu %d la: %f\n",lap,b[tam][j]/b[tam][j-1]);
		printf("vector rieng tuong ung la: ");
		for (i=1;i<=n;i++)
			printf("%f ",b[i][j]/temp);
		s=0;printf("\n");
		for (i=1;i<=n;i++) s+=pow(b[i][j],2);
		for (i=1;i<=n;i++)	//tinh ma tran moi de tinh tri rieng tiep theo
			for (k=1;k<=n;k++)
				a[i][k]-=b[tam][j]/b[tam][j-1]/s*b[i][j]*b[k][j];
		break;
		case 3:
		j++;
		for (i=1;i<n;i++)
			{s=0;
			for (k=1;k<=n;k++)
				s+=a[i][k]*b[k][j-1];
			b[i][j]=s;}
		j++;
		for (i=1;i<n;i++)
			{s=0;
			for (k=1;k<=n;k++)
				s+=a[i][k]*b[k][j-1];
			b[i][j]=s;}
		printf("\ngia tri rieng thu %d la: %f ",lap,sqrt(b[1][j]/b[1][j-2]));		//in vector rieng va tri rieng
		printf("vector rieng tuong ung la: ");
		for (i=1;i<=n;i++) printf("%f ",b[i][j]+sqrt(b[1][j]/b[1][j-2])*b[i][j-1]);
		printf("\ngia tri rieng thu %d la: %f ",lap+1,-sqrt(b[1][j]/b[1][j-2]));
		printf("vector rieng tuong ung la: ");
		for (i=1;i<=n;i++) printf("%f ",b[i][j]-sqrt(b[1][j]/b[1][j-2])*b[i][j-1]);
		printf("\ngap gia tri rieng troi trai dau, chuong trinh se dung tai day");break;
		case 4:
		j++;
		for (i=1;i<n;i++)
			{s=0;
			for (k=1;k<=n;k++)
				s+=a[i][k]*b[k][j-1];
			b[i][j]=s;}
		j++;
		for (i=1;i<n;i++)
			{s=0;
			for (k=1;k<=n;k++)
				s+=a[i][k]*b[k][j-1];
			b[i][j]=s;}
		x=1;		//tinh cac he so cua phuong trinh
		y=(b[1][j]*b[2][j-2]-b[1][j-2]*b[2][j])/(b[1][j-2]*b[2][j-1]-b[1][j-1]*b[2][j-2]);
		z=(b[1][j-1]*b[2][j]-b[1][j]*b[2][j-1])/(b[1][j-2]*b[2][j-1]-b[1][j-1]*b[2][j-2]);
		float delta=y*y-4*x*z;		//tinh delta
		printf("%f %f %f %f",x,y,z,delta);
		printf("\ntri rieng thu %d la: %.2f+%.2fi\n",lap,-y/2/x,sqrt(-delta)/2/x);
		printf("vector rieng tuong ung la: ");temp=b[1][j];//in gia tri rieng va vector rieng
		for (i=1;i<=n;i++) printf ("%.2f-(%.2f)i ",(b[i][j+1]-delta*b[i][j])/temp,(sqrt(-delta)*b[i][j])/temp);printf("\n");
		printf("\ntri rieng thu %d la: %.2f-%.2fi\n",lap+1,-y/2/x,sqrt(-delta)/2/x);
		printf("vector rieng tuong ung la: ");
		for (i=1;i<=n;i++) printf ("%.2f+(%.2f)i ",(b[i][j+1]-delta*b[i][j])/temp,(sqrt(-delta)*b[i][j])/temp);printf("\n");
		printf("gap tri rieng phuc, chuong trinh se dung tai day");
		break;
	}
	if (th==3||th==4) break;
}}




//-------------------------------
// mac.cpp
// tunaProducts(2022.12.23)
// MAC法、中心差分
//-------------------------------
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

using namespace std;

int main()
{
    //文字の定義
    int i,j,k,l,m;
    const int xn=51;//X方向節点数
	const int yn=51;//Y方向節点数
    double u[xn+2][yn+1];//X方向速度
    double v[xn+1][yn+2];//Y方向速度
    double p[xn+1][yn+1];//圧力
    double dx=1.0/(double)(xn-1);//メッシュ1つのX方向長さ
    double dy=1.0/(double)(yn-1);//メッシュ1つのY方向長さ
    double udiv,vdiv;//メッシュ間の速度差
    double d[xn+1][yn+1];//連続の式の為の項目
    double divv;//udiv+vdivの絶対値の総和
    double r[xn+1][yn+1];//ポアソン方程式の右辺項
    double ua,ub,va,vb;//各メッシュの中間位置での速度
    double err;//収斂計算の許容誤差
    double pres;//逐一更新する圧力値
	double uad,vad;
	double udif,vdif;
	double umid,vmid;
    
    //初期条件
	double uwall = 1.0;//最上部が右に1m/sで移動
	double dt = 0.001; //解析刻み
    double time = 20.0;//解析終了時間
	double re = 100;   //レイノルズ数
    int lm = time / dt;//解析ステップ数(dt*lm→解析終了時間)
    int km = 100;      //ポアソン方程式を解く際の収斂計算の最大ステップ数
	double C1 = 0.5 * dy * dy/(dx * dx + dy * dy);//SOR法計算用
	double C2 = 0.5 * dx * dx/(dx * dx + dy * dy);//SOR法計算用
	double C3 = 0.5 * dy * dy/(1.0 + dy * dy/(dx * dx));//SOR法計算用
    
    //出力設定
	ofstream fv,fp;
	fv.open("macVelOut.txt");//速度の出力
	fp.open("macPreOut.txt");//圧力の出力

    
	//速度と圧力の初期化
	for(i=0;i<xn+1;i++)
	{
		for(j=0;j<yn+1;j++)
		{
			p[i][j]=0.0;
		}
	}
	for(i=0;i<xn+2;i++)
	{
		for(j=0;j<yn+1;j++)
		{
			u[i][j]=0.0;
		}
	}
	for(i=0;i<xn+1;i++)
	{
		for(j=0;j<yn+2;j++)
		{
			v[i][j]=0.0;
		}
	}
	
	//MAC法解析開始
	for(l=1;l<=lm;l++)
	{
		//境界条件(左右)
		for(j=0;j<yn+1;j++)
		{
			u[1][j]=0.0;
			u[0][j]=u[2][j];
			v[0][j]=-v[1][j];
			
			u[xn][j]=0.0;
			u[xn+1][j]=u[xn-1][j];
			v[xn][j]=-v[xn-1][j];
		}
		v[0][yn+1]=-v[1][yn+1];
		v[xn][yn+1]=-v[xn-1][yn+1];
		
		//境界条件(上下)
		for(i=0;i<xn+1;i++)
		{
			v[i][1]=0.0;
			v[i][0]=v[i][2];
			u[i][0]=-u[i][1];
			
			v[i][yn]=0.0;
			v[i][yn+1]=v[i][yn-1];
			u[i][yn]=2.0*uwall-u[i][yn-1];
		}
		u[xn+1][0]=-u[xn][0];
		u[xn+1][yn]=-u[xn][yn];
		
		divv=0.0;
        
		//ポアソン方程式を解くための準備計算
		for(i=1;i<xn;i++)
		{
			for(j=1;j<yn;j++)
			{
				udiv=(u[i+1][j]-u[i][j])/dx;
				vdiv=(v[i][j+1]-v[i][j])/dy;
				d[i][j]=udiv+vdiv;
				divv+=fabs(udiv+vdiv);
				ua=(u[i][j]+u[i+1][j]+u[i+1][j+1]+u[i][j+1])/4.0;//右上の点
				ub=(u[i][j]+u[i+1][j]+u[i+1][j-1]+u[i][j-1])/4.0;//右下の点
				va=(v[i][j]+v[i][j+1]+v[i+1][j+1]+v[i+1][j])/4.0;//右上の点
				vb=(v[i][j]+v[i][j+1]+v[i-1][j+1]+v[i-1][j])/4.0;//左上の点
				r[i][j]=-udiv*udiv-2.0*(ua-ub)*(va-vb)/dx/dy-vdiv*vdiv+1.0/dt*(udiv+vdiv);
			}
		}
		
		//ポアソン方程式を解く
		for(k=1;k<=km;k++)
		{
			err=0.0;
			//圧力に関する境界条件
			for(j=0;j<yn+1;j++)
			{
				p[0][j]=p[1][j]-1.0/re*2.0*u[2][j];
				p[xn][j]=p[xn-1][j]+1.0/re*2.0*u[xn-1][j];
			}
			
			for(i=0;i<xn+1;i++)
			{
				p[i][0]=p[i][1]-1.0/re*2.0*v[i][2];
				p[i][yn]=p[i][yn-1]+1.0/re*2.0*v[i][yn-1];
			}
			
			//SORによる収斂計算//
			for(i=1;i<xn;i++)
			{
				for(j=1;j<yn;j++)
				{
					pres=C1*(p[i+1][j]+p[i-1][j])+C2*(p[i][j+1]+p[i][j-1])-C3*r[i][j]-p[i][j];
					err+=pres*pres;
					p[i][j]=pres+p[i][j];
				}
			}
			
			if(err<=0.000005) break;
		}
        
        //X方向速度の算出
		for(i=2;i<xn;i++)
		{
			for(j=1;j<yn;j++)
			{
				vmid=(v[i][j]+v[i][j+1]+v[i-1][j+1]+v[i-1][j])/4.0;
				uad=u[i][j]*(u[i+1][j]-u[i-1][j])/2.0/dx+vmid*(u[i][j+1]-u[i][j-1])/2.0/dy;
				udif=(u[i+1][j]-2.0*u[i][j]+u[i-1][j])/dx/dx+(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dy/dy;
				u[i][j]=u[i][j]+dt*(-uad-(p[i][j]-p[i-1][j])/dx+1.0/re*udif);
			}
		}
        
		
        //Y方向速度の算出
		for(i=1;i<xn;i++)
		{
			for(j=2;j<yn;j++)
			{
				umid=(u[i][j]+u[i+1][j]+u[i+1][j-1]+u[i][j-1])/4.0;
				vad=umid*(v[i+1][j]-v[i-1][j])/2.0/dx+v[i][j]*(v[i][j+1]-v[i][j-1])/2.0/dy;
				vdif=(v[i+1][j]-2.0*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dy/dy;
				v[i][j]=v[i][j]+dt*(-vad-(p[i][j]-p[i][j-1])/dy+1.0/re*vdif);
			}
		}
	}

    
    //数値の出力
	for(i=1;i<xn;i++)
	{
		for(j=1;j<yn;j++)
		{
			fv<<double((i-0.5)*dx)<<" "<<double((j-0.5)*dy)<<" "<<(u[i][j]+u[i+1][j])/2.0<<" "<<(v[i][j]+v[i][j+1])/2.0<<endl;
			fp<<double((i-0.5)*dx)<<" "<<double((j-0.5)*dy)<<" "<<p[i][j]<<endl;
		}
	}
	
	return 0;
}

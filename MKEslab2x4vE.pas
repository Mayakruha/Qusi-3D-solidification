unit MKEslab2x4vE;

interface

uses MKEobj, MatrixSmll, Math;

const
  N=4;
  NFree=2;
  NGuk=3;
  Wcr=0;{От 0 до 1}

type
  TCreepKoef = record
          Ac:Extended;
          kT:Extended;
          k:Extended;
    end;
        
  TMKEslab2x4=class(TMKE)
    {Свойства материала}
    CreepPr:TCreepKoef;
    E,mu,beta,dtau,Error,SigmaMax,EpsZ,dExtendForce,GeneralSqG:Extended;
    Sigma0,Sigma:array of Extended;
    Temp0,Temp:array of Extended;
    Stress0,Stress:array[0..NGuk]of Extended;
    Eps,EpsCr,EpsT:array[0..NGuk] of Extended;
    GlGuk,InvGuk:array[0..(NGuk+1)*(NGuk+1)-1]of Extended;
    {-----------------------------------------}
    constructor Create;
    destructor Destroy;override;
   { function Square(x,y:array of Extended):Extended;{находит площадь элемента}
    function CreepSpeedDivSigma(Value,Temp:Extended):Extended;{Скорость ползучести}
    function SigmaI(Sigma:array of Extended):Extended;{Интенсивность напряжений}
    procedure FuncForm(var z:array of Extended;ksi,nu:Extended);override;{Функции формы}
    procedure DifFormLoc(var z:array of Extended;ksi,nu:Extended);override;{функция формы}
    procedure GukMatrixFind(var M:array of Extended;NEl:word;NGp:byte);override;
    procedure GradientMatrixFind(var M:array of Extended;NEl:word;ksi,nu:Extended);override;
              {Нахождение матрицы градиентов}
    procedure RigidMatrixAndForceFind(var M:array of Extended;NEl:word;NGP:byte);override;
    procedure SigmaFind(var M:array of Extended;NEl:word;NGp:byte);{Нахождение Напряжения в точке элемента}
    function EpsCreepFind(NEl:word;NGP:byte):Extended;
    procedure ElementFind(var M,V:array of Extended;NEl:word);override;
    procedure ForceFind(var M:array of Extended;NEl:word;NGP:byte);
    procedure Solve;override;
    procedure NextStep;{подготовка к следующему шагу}
    procedure SetArray;
    procedure GlobalForce;
    procedure ReadElement(NEl:word);
  protected
    Kf:Extended;
    DifGlob:array[0..2*N-1]of Extended;
    dRM:array[0..NFree*NFree*N*N-1]of Extended;{Матрица приращения жесткости элемента}
    ElTemp0,ElTemp:array[0..N-1]of Extended;
    T0,T,SI0,SI,VolSig0,VolSig,dNz,SqEl:Extended;
    S3:array[0..NGuk-1]of Extended;
    S4:array[0..NGuk]of Extended;
    ijl:word;
  end;

implementation

constructor TMKEslab2x4.Create;
begin
  inherited Create(4,2,3);
  NGaussPoints:=4;
  PntG[0]:=-0.5;
  PntG[1]:=-0.5;
  PntG[2]:=0.5;
  PntG[3]:=-0.5;
  PntG[4]:=0.5;
  PntG[5]:=0.5;
  PntG[6]:=-0.5;
  PntG[7]:=0.5;
  WG[0]:=1;
  WG[1]:=1;
  WG[2]:=1;
  WG[3]:=1;
end;

destructor TMKEslab2x4.Destroy;
begin
  Temp0:=Nil;
  Temp:=Nil;
  Sigma0:=Nil;
  Sigma:=Nil;
  inherited Destroy;
end;

procedure TMKEslab2x4.SetArray;
begin
  SetLength(Temp0,NAll);
  SetLength(Temp,NAll);
  SetLength(Sigma0,NElement*NGaussPoints*(NGuk+1));
  SetLength(Sigma,NElement*NGaussPoints*(NGuk+1));
end;

{function TMKEslab2x4.Square(x,y:array of Extended):Extended;
begin
  Result:=(x[1]*y[2]-x[2]*y[1]+x[0]*y[1]-x[1]*y[0]+x[2]*y[3]-x[3]*y[2]+
   x[3]*y[0]-x[0]*y[3])/2;
end;  }

function TMKEslab2x4.CreepSpeedDivSigma(Value,Temp:Extended):Extended;
begin
  if Value<>0 then Result:=Power(Value/CreepPr.Ac*exp(CreepPr.kT*Temp),1/CreepPr.k)/Value
  else Result:=0;
end;

function TMKEslab2x4.SigmaI(Sigma:array of Extended):Extended;
begin
  Result:=Sqrt(((Sigma[0]-Sigma[1])*(Sigma[0]-Sigma[1])+(Sigma[0]-Sigma[3])*(Sigma[0]-Sigma[3])+
    (Sigma[1]-Sigma[3])*(Sigma[1]-Sigma[3]))/2+3*Sigma[2]*Sigma[2]);
end;

procedure TMKEslab2x4.FuncForm(var z:array of Extended;ksi,nu:Extended);
begin
  z[0]:=(1-ksi)*(1-nu)/4;
  z[1]:=(1+ksi)*(1-nu)/4;
  z[2]:=(1+ksi)*(1+nu)/4;
  z[3]:=(1-ksi)*(1+nu)/4;
end;

procedure TMKEslab2x4.DifFormLoc(var z:array of Extended;ksi,nu:Extended);
begin
  z[0]:=-(1-nu)/4;
  z[1]:=(1-nu)/4;
  z[2]:=(1+nu)/4;
  z[3]:=-(1+nu)/4;
  z[4]:=-(1-ksi)/4;
  z[5]:=-(1+ksi)/4;
  z[6]:=(1+ksi)/4;
  z[7]:=(1-ksi)/4;
end;

procedure TMKEslab2x4.GukMatrixFind(var M:array of Extended;NEl:word;NGp:byte);
{need - ElTemp0[], ElTemp[]}
var
  i,j:byte;
begin
  FuncForm(FormVector,ksi,nu);
  T0:=ScalarVector(FormVector,ElTemp0);
  T:=ScalarVector(FormVector,ElTemp);
  for i:=0 to NGuk do
    begin
      Stress[i]:=Sigma[NEl*NGaussPoints*(NGuk+1)+NGp*(NGuk+1)+i];
      Stress0[i]:=Sigma0[NEl*NGaussPoints*(NGuk+1)+NGp*(NGuk+1)+i];
    end;
  SI:=SigmaI(Stress);
  SI0:=SigmaI(Stress0);
  if SI=SI0 then SI:=SI0+0.5*beta/E*abs(T-T0);

  VolSig0:=(Stress0[0]+Stress0[1]+Stress0[3])/3;
  Stress[0]:=Stress0[0]-VolSig0;
  Stress[1]:=Stress0[1]-VolSig0;
  Stress[2]:=Stress0[2];
  Stress[3]:=Stress0[3]-VolSig0;

  GlGuk[0]:=1/E+dtau*(1-Wcr)*CreepSpeedDivSigma(SI,T0);
  GlGuk[1]:=-mu/E-dtau*(1-Wcr)*CreepSpeedDivSigma(SI,T0)/2;
  GlGuk[2]:=0;
  GlGuk[10]:=2*(1+mu)/E+3*dtau*(1-Wcr)*CreepSpeedDivSigma(SI,T0);

  GlGuk[3]:=GlGuk[1];
  GlGuk[5]:=GlGuk[0];
  GlGuk[6]:=0;
  GlGuk[7]:=GlGuk[1];
  GlGuk[11]:=0;
  GlGuk[15]:=GlGuk[0];

  if SI0>0 then
    begin
      if SI0<>SI then Kf:=dtau*(1-Wcr)*9/4/SI0*(CreepSpeedDivSigma(SI,T0)-CreepSpeedDivSigma(SI0,T0))/(SI-SI0)
      else Kf:=dtau*(1-Wcr)*9/4/SI0*(1/CreepPr.k-1)*CreepSpeedDivSigma(SI0,T0)/SI0;

      GlGuk[0]:=GlGuk[0]+Kf*Stress[0]*Stress[0];
      GlGuk[1]:=GlGuk[1]+Kf*Stress[0]*Stress[1];
      GlGuk[2]:=2*Kf*Stress[0]*Stress[2];
      GlGuk[3]:=GlGuk[3]+Kf*Stress[0]*Stress[3];

      GlGuk[5]:=GlGuk[5]+Kf*Stress[1]*Stress[1];
      GlGuk[6]:=2*Kf*Stress[1]*Stress[2];
      GlGuk[7]:=GlGuk[7]+Kf*Stress[1]*Stress[3];

      GlGuk[10]:=GlGuk[10]+4*Kf*Stress[2]*Stress[2];
      GlGuk[11]:=2*Kf*Stress[2]*Stress[3];

      GlGuk[15]:=GlGuk[15]+Kf*Stress[3]*Stress[3];
    end;
  GlGuk[4]:=GlGuk[1];
  GlGuk[8]:=GlGuk[2];
  GlGuk[9]:=GlGuk[6];
  GlGuk[12]:=GlGuk[3];
  GlGuk[13]:=GlGuk[7];
  GlGuk[14]:=GlGuk[11];

  OpositeMatrix(GlGuk,InvGuk,NGuk+1);
  for i:=0 to NGuk-1 do
    for j:=0 to NGuk-1 do M[i*NGuk+j]:=InvGuk[i*(NGuk+1)+j];
end;

procedure TMKEslab2x4.GradientMatrixFind(var M:array of Extended;  {need xp[],yp[]}
          NEl:word;ksi,nu:Extended);
var
  i,j,j1:byte;
begin
  DifFormGlob(DifGlob,xp,yp,ksi,nu);
  for j:=0 to N-1 do
    for j1:=0 to NFree-1 do
      for i:=0 to NGuk-1 do
        begin
          if (j1=0)and(i=0) then Varabl:=Element(0,j,DifGlob,2,N)
          else if (j1=1)and(i=(NGuk-1)) then Varabl:=Element(0,j,DifGlob,2,N)
          else if (j1=1)and(i=1) then Varabl:=Element(1,j,DifGlob,2,N)
          else if (j1=0)and(i=(NGuk-1)) then Varabl:=Element(1,j,DifGlob,2,N)
          else Varabl:=0;
          SetElement(Varabl,i,j*NFree+j1,M,NGuk,N*NFree);
        end;
end;

procedure TMKEslab2x4.RigidMatrixAndForceFind(var M:array of Extended;NEl:word;NGP:byte);
var
  i,j:byte;
begin
  GukMatrixFind(Guk,NEl,NGP); {found - FormVector[], T0, T, SI0, Stress[], InvGuk[]}
  Eps[0]:=beta*(T-T0)+dtau*3/2*CreepSpeedDivSigma(SI0,T0)*Stress[0];
  Eps[1]:=beta*(T-T0)+dtau*3/2*CreepSpeedDivSigma(SI0,T0)*Stress[1];
  Eps[2]:=dtau*3*CreepSpeedDivSigma(SI0,T0)*Stress[2];
  Eps[3]:=beta*(T-T0)+dtau*3/2*CreepSpeedDivSigma(SI0,T0)*Stress[3];

  if Global_B=Nil then GradientMatrixFind(Matrix_B,NEl,PntG[2*NGP],PntG[2*NGP+1])
  else
    for i:=0 to NGuk-1 do
      for j:=0 to N*NFree-1 do Matrix_B[i*N*NFree+j]:=
       Global_B[NEl*NGaussPoints*NGuk*N*NFree+NGP*NGuk*N*NFree+i*N*NFree+j];
  Transpon(Matrix_B,NGuk,NFree*N,Matrix_BT);
  MatrixXMatrix(Matrix_BT,Guk,NFree*N,NGuk,NGuk,Matrix_BE);
  MatrixXMatrix(Matrix_BE,Matrix_B,NFree*N,NGuk,NFree*N,dRM);

  MatrixXMatrix(Guk,Eps,NGuk,NGuk,1,S3);
  for i:=0 to NGuk-1 do S3[i]:=S3[i]-InvGuk[i*(NGuk+1)+NGuk]*(EpsZ-Eps[NGuk])+0*Stress0[i]; {невязка}
  MatrixXMatrix(Matrix_BT,S3,NFree*N,NGuk,1,VectorEl);

  for i:=0 to NFree*N-1 do
    begin
      for j:=0 to NFree*N-1 do M[i*(NFree*N+1)+j]:=dRM[i*NFree*N+j];
      M[i*(NFree*N+1)+NFree*N]:=VectorEl[i];
    end;
end;

procedure TMKEslab2x4.ForceFind(var M:array of Extended;NEl:word;NGP:byte);
var
  i,l:byte;
begin
  GukMatrixFind(Guk,NEl,NGP); {found - FormVector[], T0,T,SI0, Stress[], InvGuk[]}
  Eps[0]:=beta*(T-T0)+dtau*3/2*CreepSpeedDivSigma(SI0,T0)*Stress[0];
  Eps[1]:=beta*(T-T0)+dtau*3/2*CreepSpeedDivSigma(SI0,T0)*Stress[1];
  Eps[2]:=dtau*3*CreepSpeedDivSigma(SI0,T0)*Stress[2];
  Eps[3]:=beta*(T-T0)+dtau*3/2*CreepSpeedDivSigma(SI0,T0)*Stress[3];

  if Global_B=Nil then GradientMatrixFind(Matrix_B,NEl,PntG[2*NGP],PntG[2*NGP+1])
   else
     for i:=0 to NGuk-1 do
       for l:=0 to N*NFree-1 do Matrix_B[i*N*NFree+l]:=
         Global_B[NEl*NGaussPoints*NGuk*N*NFree+NGP*NGuk*N*NFree+i*N*NFree+l];
  Transpon(Matrix_B,NGuk,NFree*N,Matrix_BT);
  MatrixXMatrix(Guk,Eps,NGuk,NGuk,1,S3);
  for i:=0 to NGuk-1 do S3[i]:=S3[i]-InvGuk[i*(NGuk+1)+NGuk]*(EpsZ-Eps[NGuk])+0*Stress0[i]; {невязка}
  MatrixXMatrix(Matrix_BT,S3,NFree*N,NGuk,1,M);
end;

procedure TMKEslab2x4.ReadElement(NEl:word);
var
  i:byte;
begin
  for i:=0 to N-1 do
    begin
       ElTemp0[i]:=Temp0[ElementPoint[NEl*N+i]];
       ElTemp[i]:=Temp[ElementPoint[NEl*N+i]];
    end;
end;

procedure TMKEslab2x4.ElementFind(var M,V:array of Extended;NEl:word);
var
  i:byte;
begin
  for i:=0 to N-1 do
    begin
      ElTemp0[i]:=Temp0[ElementPoint[NEl*N+i]];
      ElTemp[i]:=Temp[ElementPoint[NEl*N+i]];
    end;
  inherited ElementFind(M,V,NEl);
end;

procedure TMKEslab2x4.SigmaFind(var M:array of Extended;NEl:word;NGp:byte); {need ElTemp0[], ElTemp[], ForceEl[]}
var
  i,j:word;
begin
  GukMatrixFind(Guk,NEl,NGp); {found - FormVector[], T, T0, Stress[], InvGuk[]}
  for i:=0 to NGuk-1 do
    for j:=0 to N*NFree-1 do Matrix_B[i*N*NFree+j]:=
      Global_B[NEl*NGaussPoints*NGuk*N*NFree+NGp*NGuk*N*NFree+i*N*NFree+j];
  MatrixXMatrix(Matrix_B,ForceEl,NGuk,N*NFree,1,S3);

  for j:=0 to NGuk-2 do Eps[j]:=S3[j]-beta*(T-T0)-dtau*3/2*CreepSpeedDivSigma(SI0,T0)*Stress[j];
  Eps[NGuk-1]:=S3[NGuk-1]-dtau*3*CreepSpeedDivSigma(SI0,T0)*Stress[NGuk-1];
  Eps[NGuk]:=EpsZ-beta*(T-T0)-dtau*3/2*CreepSpeedDivSigma(SI0,T0)*Stress[NGuk];
  MatrixXMatrix(InvGuk,Eps,NGuk+1,NGuk+1,1,M);
end;

function TMKEslab2x4.EpsCreepFind(NEl:word; NGp:byte):Extended; {need ElTemp0[], ElTemp[]}
var
  i:byte;
begin
  FuncForm(FormVector,ksi,nu);
  T0:=ScalarVector(FormVector,ElTemp0);
  T:=ScalarVector(FormVector,ElTemp);
  for i:=0 to NGuk do
    Stress0[i]:=Sigma0[NEl*NGaussPoints*(NGuk+1)+NGp*(NGuk+1)+i];
  SI0:=SigmaI(Stress0);
  VolSig0:=(Stress0[0]+Stress0[1]+Stress0[3])/3;
  Stress[0]:=Stress0[0]-VolSig0;
  Stress[1]:=Stress0[1]-VolSig0;
  Stress[2]:=Stress0[2];
  Stress[3]:=Stress0[3]-VolSig0;
  for i:=0 to NGuk do
    Eps[i]:=Wcr*3/2*CreepSpeedDivSigma(SI0,T0)*Stress[i];

  for i:=0 to NGuk do
    Stress0[i]:=Sigma[NEl*NGaussPoints*(NGuk+1)+NGp*(NGuk+1)+i];
  SI0:=SigmaI(Stress0);
  VolSig0:=(Stress0[0]+Stress0[1]+Stress0[3])/3;
  Stress[0]:=Stress0[0]-VolSig0;
  Stress[1]:=Stress0[1]-VolSig0;
  Stress[2]:=Stress0[2];
  Stress[3]:=Stress0[3]-VolSig0;
  for i:=0 to NGuk do
    Eps[i]:=Eps[i]+(1-Wcr)*3/2*CreepSpeedDivSigma(SI0,T)*Stress[i];

  Result:=sqrt(2/3*(Eps[0]*Eps[0]+Eps[1]*Eps[1]+Eps[3]*Eps[3]+2*Eps[2]*Eps[2]));
end;

procedure TMKEslab2x4.Solve;
var
  i,j,l:word;
begin
  dNz:=0;
  Error:=0;
  GeneralSqG:=0;
  SigmaMax:=0;
  for i:=0 to NElement-1 do
    for j:=0 to N-1 do
      begin
        for l:=0 to NGuk do S4[l]:=Sigma[N*(NGuk+1)*i+(NGuk+1)*j+l];
        SI:=SigmaI(S4);
        if SI>SigmaMax then
          begin
            SigmaMax:=SI;
            ijl:=N*(NGuk+1)*i+(NGuk+1)*j;
          end;
      end;
  inherited Solve;
  for i:=0 to NElement-1 do
    begin
      for j:=0 to N-1 do
        begin
          xp[j]:=MassivPoint[2*ElementPoint[i*N+j]];
          yp[j]:=MassivPoint[2*ElementPoint[i*N+j]+1];
          ElTemp0[j]:=Temp0[ElementPoint[i*N+j]];
          ElTemp[j]:=Temp[ElementPoint[i*N+j]];
          for l:=0 to NFree-1 do
            ForceEl[NFree*j+l]:=Vector.Massiv[NFree*(ElementPoint[i*N+j])+l];
        end;
      for j:=0 to NGaussPoints-1 do
        begin
          SigmaFind(S4,i,j);
          for l:=0 to NGuk do Sigma[NGaussPoints*(NGuk+1)*i+(NGuk+1)*j+l]:=S4[l]+Sigma0[NGaussPoints*(NGuk+1)*i+(NGuk+1)*j+l];
          Jacobian(Jacob,PntG[2*j],PntG[2*j+1],xp,yp);
          Delta:=DefMatrix(Jacob,2);
          dNz:=dNz+S4[NGuk]*WG[j]*Delta;
          GeneralSqG:=GeneralSqG+InvGuk[NGuk*(NGuk+1)+NGuk]*WG[j]*Delta;
        end;
    end;
  for l:=0 to NGuk do S4[l]:=Sigma[ijl+l];
  SI:=SigmaI(S4);
  if SigmaMax=0 then Error:=1
  else Error:=abs(1-SI/SigmaMax);
  EpsZ:=EpsZ+(dExtendForce-dNz)/GeneralSqG;
end;

procedure TMKEslab2x4.GlobalForce;
var
  i,j,l:word;
begin
  for i:=0 to NFree*NAll-1 do Force.Massiv[i]:=0;
  for j:=0 to NElement-1 do
    begin
      for i:=0 to N-1 do
        begin
          xp[i]:=MassivPoint[2*ElementPoint[N*j+i]];
          yp[i]:=MassivPoint[2*ElementPoint[N*j+i]+1];
          ElTemp0[i]:=Temp0[ElementPoint[j*N+i]];
          ElTemp[i]:=Temp[ElementPoint[j*N+i]];
        end;
      GaussIntegSquare(ForceEl,j,ForceFind);
      for i:=0 to N-1 do
        for l:=0 to NFree-1 do
          begin
            NPoint:=ElementPoint[j*N+i];
            NumPoint[i*NFree+l]:=NFree*NPoint+l;
            Force.Massiv[NFree*NPoint+l]:=Force.Massiv[NFree*NPoint+l]+ForceEl[NFree*i+l]
          end;
    end;
end;

procedure TMKEslab2x4.NextStep;
var
  i,j,l:word;
begin
  for i:=0 to NAll-1 do Temp0[i]:=Temp[i];
  for i:=0 to NElement-1 do
    for j:=0 to N-1 do
      for l:=0 to NGuk do Sigma0[N*(NGuk+1)*i+(NGuk+1)*j+l]:=Sigma[N*(NGuk+1)*i+(NGuk+1)*j+l];
end;

end.

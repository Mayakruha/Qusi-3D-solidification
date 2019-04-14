unit MainForm;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls, MKEobj, MKEslab2x4vE, Math, Matchad, MKEQ2x4;

type
  TFormCheck = class(TForm)
    Memo: TMemo;
    ButtonExit: TButton;
    Button1: TButton;
    procedure ButtonExitClick(Sender: TObject);
    procedure Button1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

  TTempMKE =class(TMKEQ2x4)
    dS:Extended;
    procedure Load(Value:Extended);
  end;

const
  BgnCalc=300{975};
  v=1.1{m/min};
  T0=1829{K};
  Tsol=1764{K};
  Nx=50{100};
  Ny=50;
  NumTempPn={1}Nx;
  Szv=7{MPa};
  dtau0=0.02{sek};
  BgnF=0.32{m};
  h0=25{mm};
  kf1=0.005;
  kf2=0.1;
  NumStep=1;

var
  FormCheck: TFormCheck;
  GlobalMKE: TMKE;
  RigidMKE: TMKEslab2x4;
  TempMKE: TTempMKE;
  {GidroMKE: TMKEgidro2x4; }
  Flag:boolean;
  UZero:array[0..Ny]of word;
  VZero:array[0..Nx]of word;
  Temp0,Temp,Sigma,Z:array of Extended;
  NRgdPoint,NRgdElem,indexF,NumWidth,MaxNumWidth,MaxNum,MinNum:word;
  DefI,DefMax,tau,Lz,Variable,dtau,dtauMin,Si,VT,xf:Extended;
  NumGlobal,NumRigid,NumGlobalElem:array of word;
  FT:File of Extended;
  F,F1:TextFile;

implementation

{$R *.DFM}

procedure TFormCheck.ButtonExitClick(Sender: TObject);
begin
  FormCheck.Close;
end;

procedure Faska(s:Extended);
begin
  if s>BgnF then xf:=(s-BgnF)*h0/0.5
  else xf:=0;
  if TempMKE.MassivPoint[2*(indexF+1)]<=xf then
    begin
      if indexF>0 then TempMKE.LoadElements[indexF*(N+1)+1]:=0;
      TempMKE.LoadElements[indexF*(N+1)+2]:=0;
      indexF:=indexF+1;
    end;
end;

function Q(s:Extended):Extended;{Vt/mm^2}
begin
  Result:=-(2.5*exp(-3.5*s/v)+0.50);
end;

procedure TTempMKE.Load(Value:Extended);
var
  i:Longword;
  j:word;
 dS:Extended;
begin
  for i:=0 to NLoadElements-1 do
    begin
      for j:=0 to N-2 do if (LoadElements[i*(N+1)+1+j]>0)and(LoadElements[i*(N+1)+1+j+1]>0) then
        begin
          NPoint:=ElementPoint[LoadElements[i*(N+1)]*N+j];
          NPoint2:=ElementPoint[LoadElements[i*(N+1)]*N+j+1];
          dS:=Sqrt(Sqr(MassivPoint[2*NPoint+0]-MassivPoint[2*NPoint2+0])+
                 Sqr(MassivPoint[2*NPoint+1]-MassivPoint[2*NPoint2+1]));
          Force.Massiv[NPoint]:=Force.Massiv[NPoint]+Q(Value)*dS/2/SqArndPnt[NPoint]*dtau
            *(Vector.Massiv[NPoint]-300)/(Vector.Massiv[NumTempPn]-300);
          Force.Massiv[NPoint2]:=Force.Massiv[NPoint2]+Q(Value)*dS/2/SqArndPnt[NPoint2]*dtau
            *(Vector.Massiv[NPoint2]-300)/(Vector.Massiv[NumTempPn]-300);
        end;
      if (LoadElements[i*(N+1)+1+0]>0)and(LoadElements[i*(N+1)+1+N-1]>0) then
        begin
          NPoint:=ElementPoint[LoadElements[i*(N+1)]*N+N-1];
          NPoint2:=ElementPoint[LoadElements[i*(N+1)]*N+0];
          dS:=Sqrt(Sqr(MassivPoint[2*NPoint+0]-MassivPoint[2*NPoint2+0])+
                 Sqr(MassivPoint[2*NPoint+1]-MassivPoint[2*NPoint2+1]));
          Force.Massiv[NPoint]:=Force.Massiv[NPoint]+Q(Value)*dS/2/SqArndPnt[NPoint]*dtau
            *(Vector.Massiv[NPoint]-300)/(Vector.Massiv[NumTempPn]-300);
          Force.Massiv[NPoint2]:=Force.Massiv[NPoint2]+Q(Value)*dS/2/SqArndPnt[NPoint2]*dtau
           *(Vector.Massiv[NPoint2]-300)/(Vector.Massiv[NumTempPn]-300);
        end;
    end;
end;

procedure TFormCheck.Button1Click(Sender: TObject);
var
  i,j,l:word;
begin
  Memo.Lines.Clear;
  GlobalMKE:=TMKE.Create(4,2,4);
  GlobalMKE.NameFileDate:='Data0.dat'{'Data8_2.txt'};
  SetLength(Temp,GlobalMKE.NAll);
  SetLength(Temp0,GlobalMKE.NAll);
  SetLength(NumRigid,GlobalMKE.NAll);
  SetLength(Sigma,GlobalMKE.NElement*GlobalMKE.N*GlobalMKE.NGuk);
  SetLength(Z,GlobalMKE.NElement*GlobalMKE.N);
  for i:=0 to High(Sigma) do Sigma[i]:=0;
  {GlobalMKE.NameFileConcen:='concen8_2.txt';
  GlobalMKE.NameFilePath:='Path8_2.txt'; }
  {--------------------------Nx,Ny--------------------}
  GlobalMKE.NLoadPoints:=Nx+Ny;
  for i:=0 to Nx-1 do
    begin
      GLobalMKE.LoadPoints[i].Point:=Ny*(Nx+1)+i;
      GlobalMKE.LoadPoints[i].LType:=2;
      GlobalMKE.LoadPoints[i].Value[1]:=0;
    end;
  for i:=0 to Ny do
    begin
      GLobalMKE.LoadPoints[i].Point:=i*(Nx+1)+Nx;
      GlobalMKE.LoadPoints[i].LType:=1;
      GlobalMKE.LoadPoints[i].Value[0]:=0;
    end;
  GlobalMKE.NPPath:=Nx+Ny+1;
  for i:=0 to Nx do GlobalMKE.Path[i]:=i;
  for i:=1 to Ny do GlobalMKE.Path[Nx+i]:=(Nx+1)*i;
  {--------------------------end:Nx,Ny-----------------}

  RigidMKE:=TMKEslab2x4.Create;
  {RigidMKE.CellVolume:=4;}
  RigidMKE.MaxNAll:=(Nx+Ny)*(25+2);
  RigidMKE.MaxWidth:=(Nx+Ny)*(25+2){(Nx+2)};
  RigidMKE.MaxNElement:=(Nx+Ny)*25;
  RigidMKE.SetArray;
  RigidMKE.E:=100000;
  RigidMKE.mu:=0.3;
  RigidMKE.beta:=1.6/100000;
  RigidMKE.CreepPr.Ac:=24000;
  RigidMKE.CreepPr.kT:=0.004;
  RigidMKE.CreepPr.k:=0.2;
  SetLength(NumGlobal,RigidMKE.MaxNAll);
  SetLength(NumGlobalElem,RigidMKE.MaxNElement);

  TempMKE:=TTempMKE.Create;
  TempMKE.NameFileDate:='Data0.dat'{'Data8_2.txt'};
  TempMKE.lamda:=29{Vt/m*K}/1000;
  TempMKE.Tsol:=Tsol{K};
  TempMKE.Tlik:=1799{K};
  TempMKE.Qkr:=272{kJ/kg}*1000;
  TempMKE.ro:=7200{kg/m3}/1000000000;
  TempMKE.Cl:=500{J/kg*K};
  TempMKE.Cr:=680{J/kg*K};
  with TempMKE do Hl:=ro*((Cl+Cr)*(Tlik-Tsol)/2+Qkr){J/m3};
  TempMKE.Epsilon:=0.0001;
  TempMKE.FormatGlobalMatrix;
  {TempMKE.NameFileDistrib:='Distrib8_2.txt';}
  TempMKE.dtau:=dtau0;
{-------------------------Nx,Ny----------------------------}
  TempMKE.NLoadElements:=Nx+Ny-1;
  TempMKE.LoadElements[0]:=0;
  TempMKE.LoadElements[0+1]:=1;
  TempMKE.LoadElements[0+2]:=1;
  TempMKE.LoadElements[0+4]:=1;
  for i:=1 to Nx-1 do
    begin
      TempMKE.LoadElements[(N+1)+(i-1)*(N+1)]:=i;
      TempMKE.LoadElements[(N+1)+(i-1)*(N+1)+1]:=1;
      TempMKE.LoadElements[(N+1)+(i-1)*(N+1)+2]:=1;
    end;
  for i:=1 to Ny-1 do
    begin
      TempMKE.LoadElements[Nx*(N+1)+(i-1)*(N+1)]:=i*Nx;
      TempMKE.LoadElements[Nx*(N+1)+(i-1)*(N+1)+1]:=1;
      TempMKE.LoadElements[Nx*(N+1)+(i-1)*(N+1)+4]:=1;
    end;
{------------------------end Nx,Ny ------------------------}
  for i:=0 to TempMKE.NAll-1 do
    begin
      Temp0[i]:=T0;
      TempMKE.Vector.Massiv[i]:=T0;
      TempMKE.Force.Massiv[i]:=TempMKE.FuncQFromTemp(T0);
    end;
  AssignFile(F,'Graphic.txt');
  Rewrite(F);
  tau:=0;
  Lz:=0;
  if FileExists('prom.ext') then
    begin
      AssignFile(FT,'prom.ext');
      Reset(FT);
      Read(FT,tau);
      Lz:=tau*v/60;
      for i:=0 to GlobalMKE.NAll-1 do
        begin
          Read(FT,Temp0[i]);
          TempMKE.Vector.Massiv[i]:=Temp0[i];
          TempMKE.Force.Massiv[i]:=TempMKE.FuncQFromTemp(Temp0[i]);
          for j:=0 to NFree-1 do
            begin
              Read(FT,Variable);
              GlobalMKE.Vector.Massiv[NFree*i+j]:=Variable;
            end;
        end;
      for i:=0 to GlobalMKE.NElement-1 do
        for j:=0 to N-1 do
          begin
            for l:=0 to NGuk do Read(FT,Sigma[i*N*(NGuk+1)+j*(NGuk+1)+l]);
            Read(FT,Z[i*N+j]);
          end;
      CloseFile(FT);
    end;
  repeat
    NRgdPoint:=0;
    for i:=0 to NumStep-1 do
      begin
        {Faska(Lz+i*dtau0*v/60);   }
        TempMKE.Load(Lz+i*dtau0*v/60);
        TempMKE.Solve;
      end;
    for i:=0 to GlobalMKE.NAll-1 do
      begin
        Temp[i]:=TempMKE.Vector.Massiv[i];
        if Temp0[i]<=Tsol then NRgdPoint:=NRgdPoint+1;
      end;
    if NRgdPoint>BgnCalc then
      begin
        {---------------------------forming rigid zone --------------------------}
        RigidMKE.dtau:=dtau0*NumStep;
        RigidMKE.NAll:=NRgdPoint;
        NRgdPoint:=0;
        for i:=0 to GlobalMKE.NAll-1 do if Temp0[i]<=Tsol then
          begin
            NumGlobal[NRgdPoint]:=i;
            NumRigid[i]:=NRgdPoint;
            for j:=0 to 1 do RigidMKE.MassivPoint[2*NRgdPoint+j]:=GlobalMKE.MassivPoint[2*i+j];
            RigidMKE.Temp0[NRgdPoint]:=Temp0[i];
            RigidMKE.Temp[NRgdPoint]:=Temp[i];
            NRgdPoint:=NRgdPoint+1
          end;
        NRgdElem:=0;
        MaxNumWidth:=0;
        for i:=0 to GlobalMKE.NElement-1 do
          begin
            Flag:=True;
            for j:=0 to N-1 do if Temp0[GlobalMKE.ElementPoint[i*N+j]]>Tsol then Flag:=False;
            if Flag then
              begin
                NumGlobalElem[NRgdElem]:=i;
                for j:=0 to N-1 do
                  begin
                    RigidMKE.ElementPoint[NRgdElem*N+j]:=NumRigid[GlobalMKE.ElementPoint[i*N+j]];
                    if (j=0)or(MinNum>RigidMKE.ElementPoint[NRgdElem*N+j]) then
                       MinNum:=RigidMKE.ElementPoint[NRgdElem*N+j];
                    if (j=0)or(MaxNum<RigidMKE.ElementPoint[NRgdElem*N+j]) then
                       MaxNum:=RigidMKE.ElementPoint[NRgdElem*N+j];
                    for l:=0 to NGuk do
                      begin
                        RigidMKE.Sigma0[NRgdElem*N*(NGuk+1)+j*(NGuk+1)+l]:=Sigma[i*N*(NGuk+1)+j*(NGuk+1)+l];
                        RigidMKE.Sigma[NRgdElem*N*(NGuk+1)+j*(NGuk+1)+l]:=Sigma[i*N*(NGuk+1)+j*(NGuk+1)+l];
                      end;
                  end;
                if MaxNumWidth<(MaxNum-MinNum) then MaxNumWidth:=MaxNum-MinNum;
                NRgdElem:=NRgdElem+1;
              end;
          end;
        RigidMKE.NElement:=NRgdElem;
        RigidMKE.Width:=MaxNumWidth;
         {---------------------------solve rigid zone----------------------------}
        RigidMKE.CalcGlobal_B;
        repeat
          RigidMKE.FormatGlobalMatrix;
          for i:=0 to GlobalMKE.NLoadPoints-1 do
            if GlobalMKE.LoadPoints[i].Point=NumGlobal[NumRigid[GlobalMKE.LoadPoints[i].Point]] then
              begin
                if GlobalMKE.LoadPoints[i].LType=1 then
                  RigidMKE.ApplyOnPoint(NumRigid[GlobalMKE.LoadPoints[i].Point],0,1,0);
                if GlobalMKE.LoadPoints[i].LType=2 then
                  RigidMKE.ApplyOnPoint(NumRigid[GlobalMKE.LoadPoints[i].Point],1,1,0);
              end;
          RigidMKE.Solve;
    {      RigidMKE.GlobalForce; }
        until RigidMKE.Error<0.01;

        {----------------------------out result----------------------------------}
        for i:=0 to NRgdPoint-1 do
          for j:=0 to NFree-1 do
            GlobalMKE.Vector.Massiv[NumGlobal[i]*NFree+j]:=
            GlobalMKE.Vector.Massiv[NumGlobal[i]*NFree+j]+RigidMKE.Vector.Massiv[i*NFree+j];
        DefMax:=0;
        for i:=0 to NRgdElem-1 do
          begin
            RigidMKE.ReadElement(i);
            for j:=0 to N-1 do
              begin
                Z[NumGlobalElem[i]*N+j]:=Z[NumGlobalElem[i]*N+j]+
                   RigidMKE.EpsCreepFind(i,j);
                if DefMax<Z[NumGlobalElem[i]*N+j] then DefMax:=Z[NumGlobalElem[i]*N+j];
                for l:=0 to NGuk do Sigma[NumGlobalElem[i]*N*(NGuk+1)+j*(NGuk+1)+l]:=RigidMKE.Sigma[i*N*(NGuk+1)+j*(NGuk+1)+l];
              end;
            end;
        Writeln(F,ToFormatMatchad(Lz)+','+ToFormatMatchad(DefMax));
        GlobalMKE.OutDate(Z,'Deform','Result.mke');

        Memo.Lines.Add(FloatToStr(Lz)+' ; '+IntToStr(NRgdElem)+' ; '+FloatToStr(GlobalMKE.Vector.Massiv[0])+' ; '
        +FloatToStr(DefMax)+' ; '+FloatToStr(RigidMKE.SigmaMax)+' ; '+FloatToStr(Temp[0]));
        Application.ProcessMessages;
      end;
    tau:=tau+dtau0*NumStep;
    Lz:=v*tau/60;{m}
 {------------------------------------------------------------------------}
    AssignFile(FT,'prom.ext');
    Rewrite(FT);
    Write(FT,tau);
    for i:=0 to GlobalMKE.NAll-1 do
      begin
        Write(FT,Temp[i]);
        for j:=0 to NFree-1 do
          begin
            Variable:=GlobalMKE.Vector.Massiv[NFree*i+j];
            Write(FT,Variable);
          end;
      end;
    for i:=0 to GlobalMKE.NElement-1 do
      for j:=0 to N-1 do
        begin
          for l:=0 to NGuk-1 do Write(FT,Sigma[i*N*NGuk+j*Nguk+l]);
          Write(FT,Z[i*N+j]);
        end;
    CloseFile(FT);
    for i:=0 to GlobalMKE.NAll-1 do Temp0[i]:=Temp[i];
  until Lz>=0.82;
  {--------вывод температуры---------------------}
 { for i:=0 to GlobalMKE.NElement-1 do
    for j:=0 to GlobalMKE.N-1 do Z[i*N+j]:=Temp[GlobalMKE.ElementPoint[i*N+j]];
  {---------end:вывод температуры---------------}
  GlobalMKE.OutDate(Z,'Deform','Result.mke');

 { if GlobalMKE.Path<>Nil then
    begin
      AssignFile(F1,'Eps.txt');
      Rewrite(F1);
      for i:=0 to GlobalMKE.NPPath-1 do
        Writeln(F1,ToFormatMatchad(GlobalMKE.MassivPoint[2*GlobalMKE.Path[i]])+
         ','+ToFormatMatchad(GlobalMKE.MassivPoint[2*GlobalMKE.Path[i]+1])+
         ','+ToFormatMatchad(Z[GlobalMKE.Path[i]]));
      CloseFile(F1);
    end;      }

  CloseFile(F);
  NumGlobal:=Nil;
  NumGlobalElem:=Nil;
  RigidMKE.Destroy;
  NumRigid:=Nil;
  Temp0:=Nil;
  Temp:=Nil;
  Sigma:=Nil;
  Z:=Nil;
  GlobalMKE.Destroy;
end;

end.

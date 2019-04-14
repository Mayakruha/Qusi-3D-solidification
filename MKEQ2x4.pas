unit MKEQ2x4;

interface

uses MKEobj, MatrixSmll, Geometry;

const
  N=4;
  NFree=1;
  NGuk=2;
  Line:array[0..3,0..3]of array[0..3]of Extended=
  (((-1,-1,0,-1),(0,-1,0,0),(0,0,-1,0),(-1,0,-1,-1)),
  ((1,-1,1,0),(1,0,0,0),(0,0,0,-1),(0,-1,1,-1)),
  ((1,1,0,1),(0,1,0,0),(0,0,1,0),(1,0,1,1)),
  ((-1,1,-1,0),(-1,0,0,0),(0,0,0,1),(0,1,-1,1)));

type
  TMKEQ2x4=class(TMKE)
    lamda,Tsol,Tlik,Cr,Cl,Qkr,ro,Hl,dtau,Epsilon,dTmax:Extended;{Свойства материала}
    SqArndPnt:array of Extended;
    constructor Create;
    destructor Destroy;override;
    function Cef(Temp:Extended):Extended;
    function FuncQFromTemp(Value:Extended):Extended;
    function TempFromFuncQ(Value:Extended):Extended;
    procedure FuncForm(var z:array of Extended;ksi,nu:Extended);override;{Функции формы}
    procedure DifFormLoc(var z:array of Extended;ksi,nu:Extended);override;{функция формы}
    procedure GukMatrixFind(var M:array of Extended;NEl:word;NGp:byte);override;
    procedure GradientMatrixFind(var M:array of Extended;NEl:word;ksi,nu:Extended);override;
              {Нахождение матрицы градиентов}
    procedure ElementFind(var M,V:array of Extended;NEl:word);override;
              {Нахождение матрицы жескости и вектора нагрузки элемента}
    procedure FindQ(var M:array of Extended;NEl:word;NGp:byte);
    procedure FormatGlobalMatrix;override;
    procedure Solve;override;

  protected
    Tk,Hk,dH,eps,Qds,T0:Extended;
    XpartEl,YpartEl:array[0..N-1]of Extended;
    VectQ:array[0..NGuk-1]of Extended;
  end;

implementation

constructor TMKEQ2x4.Create;
var
  i:byte;
begin
  inherited Create(4,1,2);
  NGaussPoints:=2;
  for i:=0 to 1 do WG[i]:=0.5;
end;

destructor TMKEQ2x4.Destroy;
begin
  if SqArndPnt<>Nil then SqArndPnt:=Nil;
  inherited Destroy;
end;

function TMKEQ2x4.Cef(Temp:Extended):Extended;
begin
  if Temp<Tsol then Result:=Cr
  else if Temp>Tlik then Result:=Cl
  else Result:=(Qkr+Cr*(Tlik-Temp)+Cl*(Temp-Tsol))/(Tlik-Tsol);
end;

function TMKEQ2x4.FuncQFromTemp(Value:Extended):Extended;
begin
  if Value<Tsol then Result:=(Value-Tsol)*ro*Cr
  else if Value>Tlik then Result:=(Value-Tlik)*ro*Cl+Hl
  else Result:=ro*(Value-Tsol)/(Tlik-Tsol)*(Qkr+(Tlik-Tsol)*Cr+
      (Cl-Cr)*(Value-Tsol)/2);
end;

function TMKEQ2x4.TempFromFuncQ(Value:Extended):Extended;
begin
  if Value<0 then Result:=Value/ro/Cr+Tsol
  else if Value>Hl then Result:=(Value-Hl)/ro/Cl+Tlik
  else
    begin
      Tk:=Value/Hl*(Tlik-Tsol)+Tsol;
      repeat
        Hk:=FuncQFromTemp(Tk);
        dH:=Value-Hk;
        eps:=abs(dH/Hl);
        Tk:=dH/ro/Cef(Tk)+Tk;
      until eps<Epsilon;
      Result:=Tk;
    end;
end;

procedure TMKEQ2x4.FuncForm(var z:array of Extended;ksi,nu:Extended);
begin
  z[0]:=(1-ksi)*(1-nu)/4;
  z[1]:=(1+ksi)*(1-nu)/4;
  z[2]:=(1+ksi)*(1+nu)/4;
  z[3]:=(1-ksi)*(1+nu)/4;
end;

procedure TMKEQ2x4.DifFormLoc(var z:array of Extended;ksi,nu:Extended);
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

procedure TMKEQ2x4.GukMatrixFind(var M:array of Extended;NEl:word;NGp:byte);
begin
  M[0]:=-lamda;
  M[1]:=0;
  M[2]:=0;
  M[3]:=-lamda;
end;

procedure TMKEQ2x4.GradientMatrixFind(var M:array of Extended;
          NEl:word;ksi,nu:Extended);{need xp[], yp[]}
begin
  DifFormGlob(M,xp,yp,ksi,nu);
end;

procedure TMKEQ2x4.ElementFind(var M,V:array of Extended;NEl:word);
begin
{-----------------------------------}
end;

procedure TMKEQ2x4.FindQ(var M:array of Extended;NEl:word;NGp:byte);
{need ForceEl[]}
begin
  ksi:=PntG[2*NGp];
  nu:=PntG[2*NGp+1];
  GradientMatrixFind(Matrix_B,NEl,ksi,nu);
  GukMatrixFind(Guk,NEl,NGp);
  MatrixXMatrix(Guk,Matrix_B,NGuk,NGuk,N,Matrix_BE);
  MatrixXMatrix(Matrix_BE,ForceEl,NGuk,N,1,M);
end;

procedure TMKEQ2x4.FormatGlobalMatrix;
var
  i,j,s:integer;
begin
  SetLength(SqArndPnt,NAll);
  for i:=0 to NAll-1 do SqArndPnt[i]:=0;
  for i:=0 to NElement-1 do
    begin
      for j:=0 to N-1 do
        begin
          xp[j]:=MassivPoint[2*ElementPoint[N*i+j]];
          yp[j]:=MassivPoint[2*ElementPoint[N*i+j]+1];
        end;
      for j:=0 to N-1 do
        begin
          for s:=0 to N-1 do
            begin
              FuncForm(FormVector,Line[j,s,0],Line[j,s,1]);
              XpartEl[s]:=ScalarVector(FormVector,xp);
              YpartEl[s]:=ScalarVector(FormVector,yp);
            end;
          SqArndPnt[ElementPoint[N*i+j]]:=SqArndPnt[ElementPoint[N*i+j]]
           +SquareTriangle(XpartEl[0],YpartEl[0],XpartEl[1],YpartEl[1],XpartEl[2],YpartEl[2])
           +SquareTriangle(XpartEl[0],YpartEl[0],XpartEl[2],YpartEl[2],XpartEl[3],YpartEl[3]);
        end;
    end;
end;

procedure TMKEQ2x4.Solve;
var
  i,j,s:integer;
begin
  for i:=0 to NElement-1 do
    begin
      for j:=0 to N-1 do
        begin
          xp[j]:=MassivPoint[2*ElementPoint[N*i+j]];
          yp[j]:=MassivPoint[2*ElementPoint[N*i+j]+1];
          ForceEl[j]:=Vector.Massiv[ElementPoint[N*i+j]];
        end;
      for j:=0 to N-1 do
        begin
          for s:=0 to N-1 do
            begin
              FuncForm(FormVector,Line[j,s,0],Line[j,s,1]);
              XpartEl[s]:=ScalarVector(FormVector,xp);
              YpartEl[s]:=ScalarVector(FormVector,yp);
            end;
          for s:=1 to 2 do
            begin
              PntG[0]:=Line[j,s,0];
              PntG[1]:=Line[j,s,1];
              PntG[2]:=Line[j,s,2];
              PntG[3]:=Line[j,s,3];
              GaussIntegLine(VectQ,i,FindQ);
              Qds:=(VectQ[0]*(YpartEl[s]-YpartEl[s+1])+VectQ[1]*(XpartEl[s+1]-XpartEl[s]))
                  /sqrt((YpartEl[s]-YpartEl[s+1])*(YpartEl[s]-YpartEl[s+1])+
                 (XpartEl[s+1]-XpartEl[s])*(XpartEl[s+1]-XpartEl[s]));
              Force.Massiv[ElementPoint[N*i+j]]:=Force.Massiv[ElementPoint[N*i+j]]
                  +Qds/SqArndPnt[ElementPoint[N*i+j]]*dtau;
            end;
        end;
    end;
  dTmax:=0;
  for i:=0 to NAll-1 do
    begin
      T0:=Vector.Massiv[i];
      Vector.Massiv[i]:=TempFromFuncQ(Force.Massiv[i]);
      if abs(dTmax)<abs(Vector.Massiv[i]-T0)then dTmax:=Vector.Massiv[i]-T0;
    end;
end;

end.

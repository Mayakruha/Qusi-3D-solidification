unit MKEobj2x4;

interface

uses MKEobj, MatrixSmll;

const
  N=4;
  NFree=2;
  NGuk=3;

type
  TMKE2x4=class(TMKE)
    {Свойства материала}
    E,mu:Extended;
    {-----------------------------------------}
    constructor Create;
    destructor Destroy;override;
    function Square(x,y:array of Extended):Extended;{находит площадь элемента}
    procedure FuncForm(var z:array of Extended;ksi,nu:Extended);override;{Функции формы}
    procedure DifFormLoc(var z:array of Extended;ksi,nu:Extended);override;{функция формы}
    procedure GukMatrixFind(var M:array of Extended;NEl:word);override;
    procedure GradientMatrixFind(var M:array of Extended;NEl:word;ksi,nu:Extended);override;
              {Нахождение матрицы градиентов}
    procedure ElementFind(var M,V:array of Extended;NEl:word);override;
              {Нахождение матрицы жескости и вектора нагрузки элемента}
    procedure SigmaFind(var Sigma:array of Extended;Eps_dop:array of Extended;
                        NEl:word;ksi,nu:Extended);{Нахождение Напряжения в точке элемента}

  protected
    DifGlob:array[0..2*N-1]of Extended;
    Eps:array[0..NGuk-1] of Extended;{}
    dRM:array[0..NFree*NFree*N*N-1]of Extended;{Матрица приращения жесткости элемента}
  end;

implementation

constructor TMKE2x4.Create;
begin
  inherited Create(4,2,3);
end;

destructor TMKE2x4.Destroy;
begin
  inherited Destroy;
end;

function TMKE2x4.Square(x,y:array of Extended):Extended;
begin
  Result:=(x[1]*y[2]-x[2]*y[1]+x[0]*y[1]-x[1]*y[0]+x[2]*y[3]-x[3]*y[2]+
   x[3]*y[0]-x[0]*y[3])/2;
end;

procedure TMKE2x4.FuncForm(var z:array of Extended;ksi,nu:Extended);
begin
  z[0]:=(1-ksi)*(1-nu)/4;
  z[1]:=(1+ksi)*(1-nu)/4;
  z[2]:=(1+ksi)*(1+nu)/4;
  z[3]:=(1-ksi)*(1+nu)/4;
end;

procedure TMKE2x4.DifFormLoc(var z:array of Extended;ksi,nu:Extended);
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

procedure TMKE2x4.GukMatrixFind(var M:array of Extended;NEl:word);
begin
  M[0]:=E/(1-mu*mu);
  M[1]:=E*mu/(1-mu*mu);
  M[2]:=0;
  M[3]:=E*mu/(1-mu*mu);
  M[4]:=E/(1-mu*mu);
  M[5]:=0;
  M[6]:=0;
  M[7]:=0;
  M[8]:=E/2/(1+mu);
end;

procedure TMKE2x4.GradientMatrixFind(var M:array of Extended;
          NEl:word;ksi,nu:Extended);
var
  i,j,j1:byte;
begin
  for i:=0 to N-1 do
    begin
      xp[i]:=MassivPoint[2*ElementPoint[NEl*N+i]];
      yp[i]:=MassivPoint[2*ElementPoint[NEl*N+i]+1];
    end;
  DifFormGlob(DifGlob,xp,yp,ksi,nu);
  for j:=0 to N-1 do
    for j1:=0 to NFree-1 do
      for i:=0 to NGuk-1 do
        begin
          if ((j1=0)and(i=0))or((j1=1)and(i=2)) then Varabl:=Element(0,j,DifGlob,2,N)
          else if ((j1=0)and(i=2))or((j1=1)and(i=1)) then Varabl:=Element(1,j,DifGlob,2,N)
          else Varabl:=0;
          SetElement(Varabl,i,j*NFree+j1,M,NGuk,N*NFree);
        end;
end;

procedure TMKE2x4.ElementFind(var M,V:array of Extended;NEl:word);
var
  i:word;
begin
  GaussInteg(M,NEl,RigidMatrixFind);

  for i:=0 to NFree*N-1 do V[i]:=0;
end;

procedure TMKE2x4.SigmaFind(var Sigma:array of Extended;Eps_dop:array of Extended;
                            NEl:word;ksi,nu:Extended);
var
  i,j:word;
begin
  for i:=0 to N-1 do
    begin
       xp[i]:=MassivPoint[2*(ElementPoint[NEl*N+i])];
       yp[i]:=MassivPoint[2*(ElementPoint[NEl*N+i])+1];
       for j:=0 to NFree-1 do ForceEl[NFree*i+j]:=Vector.Massiv[NFree*(ElementPoint[NEl*N+i])+j];
    end;
  GukMatrixFind(Guk,NEl);
  GradientMatrixFind(Matrix_B,NEl,ksi,nu);
  MatrixXMatrix(Matrix_B,ForceEl,NGuk,N*NFree,1,Eps);
  for j:=0 to NGuk-1 do Eps[j]:=Eps[j]-Eps_dop[j];
  MatrixXMatrix(Guk,Eps,NGuk,NGuk,1,Sigma);
end;

end.

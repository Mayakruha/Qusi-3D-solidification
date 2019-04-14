unit Geometry;

interface

type
   set255=set of 0..255;

{Вычисляет угол между вектороми. Угол отмеряется от вектора с индексом 0
против часовой стрелки.}
Function AngleVec(x0,y0,x1,y1:Extended):Extended;

{Определяет принадлежит ли точка с координатами Х0,У0 фигуре ограниченой
линией состоящих из точек координаты кооторых записаны в массивах.
Координаты узлов нумируются против часовой стрелки.}
Function PointBelong(x0,y0:Extended;x,y:array of Extended):Boolean;

{Выдыет значение поля в точке элемента с координатами х,у}
Function ValueXY_Element(x,y:Extended; xk,yk,zk:array of Extended;n:integer):Extended;

{Вычисляет площадь треуголника с координатами X,Y}
Function SquareTriangle(x0,y0,x1,y1,x2,y2:Extended):Extended;

{Находит следующее значение во множестве после Value}
Function Next(St:set255;Value:byte):byte;

{Находит предыдущме значение во множестве перед Value}
Function Prior(St:set255;Value:byte):byte;

{Вычисляет площадь многоугольника}
Function Square(x,y:array of Extended):Extended;

{Равен ли N-ый бит Х единице}
Function Bit(x:byte;N:byte):Boolean;

implementation

uses Math, SysUtils, MatrixSmll;

Function AngleVec(x0,y0,x1,y1:Extended):Extended;
var
  fy0,fy1:Extended;
begin
  if (x0=0)and(y0=0)then fy0:=0
  else if (x0=0)and(y0>0) then fy0:=pi/2
  else if (x0=0)and(y0<0) then fy0:=-pi/2
  else fy0:=ArcTan2(y0,x0);

  if (x1=0)and(y1=0) then fy1:=0
  else if (x1=0)and(y1>0) then fy1:=pi/2
  else if (x1=0)and(y1<0) then fy1:=-pi/2
  else fy1:=ArcTan2(y1,x1);

  if fy1<fy0 then Result:=2*pi+fy1-fy0
  else Result:=fy1-fy0;
end;

Function PointBelong(x0,y0:Extended;x,y:array of Extended):Boolean;
var
  i,imin,n1,n2:integer;
  dx0,dy0,dx1,dy1,dx2,dy2,fy0,fy1,fy2,R,Rmin:Extended;
label
  m;
begin
  n1:=Low(x);
  n2:=High(x);
  if (Low(y)<>n1)or(High(y)<>n2) then
    begin
      Result:=False;
      goto m;
    end;
  Rmin:=(x0-x[n1])*(x0-x[n1])+(y0-y[n1])*(y0-y[n1]);
  imin:=n1;
  for i:=n1+1 to n2 do
    begin
      R:=(x0-x[i])*(x0-x[i])+(y0-y[i])*(y0-y[i]);
      if R<Rmin then
        begin
          Rmin:=R;
          imin:=i;
        end;
    end;
  dx0:=x0-x[imin];
  dy0:=y0-y[imin];
  if imin=n1 then
    begin
      dx1:=x[n2]-x[n1];
      dx2:=x[n1+1]-x[n1];
      dy1:=y[n2]-y[n1];
      dy2:=y[n1+1]-y[n1];
    end
  else if imin=n2 then
    begin
      dx1:=x[n2-1]-x[n2];
      dx2:=x[n1]-x[n2];
      dy1:=y[n2-1]-y[n2];
      dy2:=y[n1]-y[n2];
    end
  else
    begin
      dx1:=x[imin-1]-x[imin];
      dx2:=x[imin+1]-x[imin];
      dy1:=y[imin-1]-y[imin];
      dy2:=y[imin+1]-y[imin];
    end;
  fy0:=AngleVec(dx2,dy2,dx1,dy1);
  fy1:=AngleVec(dx0,dy0,dx1,dy1);
  fy2:=AngleVec(dx2,dy2,dx0,dy0);
  if (fy0>=fy1)and(fy0>=fy2) then Result:=True
  else Result:=False;
m:
end;

Function ValueXY_Element(x,y:Extended; xk,yk,zk:array of Extended; n:integer):Extended;
var
  a,a_op:array of Extended;
  c,prom:array of Extended;
  move:array of integer;
  Value:Extended;
  x0,y0,z0,i,j:integer;
begin
  SetLength(a,n*n);
  SetLength(a_op,n*n);
  SetLength(c,n);
  SetLength(prom,n);
  SetLength(move,n);
  x0:=Low(xk);
  y0:=Low(yk);
  z0:=Low(zk);
  for i:=0 to n-2 do c[i]:=Power(x,n-2-i)*Power(y,i);
  c[n-1]:=1;
  for i:=0 to n-1 do
    begin
      for j:=0 to n-2 do
        begin
          Value:=Power(xk[i],n-2-j)*Power(yk[i],j);
          SetElement(Value,i,j,a,n,n);
        end;
        SetElement(1,i,n-1,a,n,n);
    end;
  if ZeroInDiag(a,n) then
    begin
      ZeroRemovDiag(a,move,n,0);
      ChangeVector(c,move,0);
    end;
  OpositeMatrix(a,a_op,n);
  MatrixSquareXVector(a_op,zk,prom);
  Result:=ScalarVector(prom,c);
  a:=Nil;
  a_op:=Nil;
  c:=Nil;
  prom:=Nil;
  move:=Nil;
end;

Function SquareTriangle(x0,y0,x1,y1,x2,y2:Extended):Extended;
begin
  Result:=(x1*y2-x2*y1+x2*y0-x0*y2+x0*y1-y0*x1)/2;
end;

Function Next(St:set255; Value:byte):byte;
var
  i:byte;
begin
  Result:=Value;
  i:=Value;
  while Result=Value do
    begin
      i:=i+1;
      if i in St then Result:=i;
    end;
end;

Function Prior(St:set255; Value:byte):byte;
var
  i:byte;
begin
  Result:=Value;
  i:=Value;
  while Result=Value do
    begin
      i:=i-1;
      if i in St then Result:=i;
    end;
end;

Function Square(x,y:array of Extended):Extended;
label
  m;
var
  Xtr,Ytr:array [0..2]of Extended;
  i:array[0..2]of byte;
  indx:set of 0..255;
  Np,j:byte;
begin
  Result:=0;
  Np:=High(x)+1;
  indx:=[0..Np-1];
  repeat
    if 0 in indx then i[0]:=0
    else i[0]:=Next(indx,0);
m:  i[1]:=Next(indx,i[0]);
    i[2]:=Next(indx,i[1]);
    for j:=0 to 2 do
      begin
        Xtr[j]:=x[i[j]];
        Ytr[j]:=y[i[j]];
      end;
    if ((Xtr[2]-Xtr[1])*(Ytr[0]-Ytr[1])+(Ytr[2]-Ytr[1])*(Xtr[1]-Xtr[0]))<0 then
      begin
        i[0]:=i[1];
        goto m;
      end;
    Np:=Np-1;
    Exclude(indx,i[1]);
    Result:=Result+SquareTriangle(Xtr[0],Ytr[0],Xtr[1],Ytr[1],Xtr[2],Ytr[2]);
  until Np<3;
end;

Function Bit(x:byte;N:byte):Boolean;
var
  i,Two,Value:byte;
begin
  Two:=1;
  for i:=1 to N-1 do Two:=Two*2;
  Value:=x mod (2*Two);
  Value:=Value div Two ;
  if Value=1 then Result:=True
  else Result:=False;
end;

end.

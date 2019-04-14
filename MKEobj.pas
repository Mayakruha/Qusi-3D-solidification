unit MKEobj;

interface

uses Matrix, Matchad, MatrixSmll, SysUtils;

type
  TFindMatrix=procedure(var M:array of Extended;NEl:word;NGp:byte) of object;

  TLoadPoints=record
    Point:word;
    LType:byte;
    Value:array[0..2]of Extended;
  end;

  TMKE=class
    N:word;{число узлов в конечном элементе}
    NFree:word;{число степеней свободы узла элемента}
    NGuk:word;{размерность матрицы свойств материала}
    CellVolume:byte;{размер ячейки для вещественного числа}
    NElement:word;{число элементов}
    MassivPoint:array of Extended;{Массив координат узлов}
    ElementPoint:array of word;{Массив номеров узлов для элементов}
    ElementMaterial:array of byte;{Массив номеров материала для элементов}
    LoadElements:array of word;{Массив элементов нагруженных распределенной нагрузкой}
    LoadPoints:array of TLoadPoints;{Maccив нагруженных точек}
    Path:array of word;{Массив номеров узлов лежащих на пути}
    MatrixRigid:TMatrix;{Общая матрица жесткости}
    Vector:TMatrix;{Глобальный вектор перемещений узлов}
    Force:TMatrix;{Глобальный вектор нагрузки}
    Global_B:array of Extended; {Массив матриц разностей по гауссовым точкам элементов}
    PntG:array of Extended; {Массив координат гауссовых точек}
    WG:array of Extended; {Массив весовых коэффициенто гауссовых точек}
{-------------------------------------------------------}
    xp,yp:array of Extended;{Координаты Х,Y узлов элемента}
    Jacob,Inv_Jacob:array[0..3]of Extended;{Якобиан,Обратный Якобиан}
    DifLoc,FormVector,Guk,Matrix_B,Matrix_BT,Matrix_BE,ElemMatrix:array of Extended;
    Xmax,Xmin,Ymax,Ymin,Zmax,Zmin:Extended;

  private
    FNameFileDate,FNameFileConcen,FNameFileDistrib,FNameFilePath:string;
    FNAll,FWidth,MaxFNAll,MaxFNElement,MaxFWidth,FNLoadElements,FNLoadPoints,FNPPath:word;
    FNGaussPoints:byte;
    {------------------ Property procedure------------------}
    procedure ReadingDate(Value:string);
    procedure FormatConcenLoad(Value:string);
    procedure FormatDistribLoad(Value:string);
    procedure FormatPath(Value:string);
    procedure SetNAll(Value:word);
    procedure SetWidth(Value:word);
    procedure SetMaxNAll(Value:word);
    procedure SetMaxNElement(Value:word);
    procedure SetNLoadElements(Value:word);
    procedure SetNLoadPoints(Value:word);
    procedure SetNPPath(Value:word);
    procedure SetMemory(Value:word);
    procedure SetNGaussPoints(Value:byte);
    procedure ReadMyFormatDate;
    procedure ReadAnsysDate;

  published
    {---------------------property----------------------------}
    {Имя файла с данными}
    property NameFileDate:string read FNameFileDate write ReadingDate;
    {Имя файла с точечным нагружением}
    property NameFileConcen:string read FNameFileConcen write FormatConcenLoad;
    {Имя файла с распределенным нагружением}
    property NameFileDistrib:string read FNameFileDistrib write FormatDistribLoad;
    {Имя файла с информацие о пути}
    property NameFilePath:string read FNameFilePath write FormatPath;
    {число узлов в расчитываемом объекте}
    property NAll:word read FNAll write SetNAll;
    {максимальная разница между номерами узлов}
    property Width:word read FWidth write SetWidth;
    {максимально возможное число узлов}
    property MaxNAll:word read MaxFNAll write SetMaxNAll;
    {максимально возможное число элементов}
    property MaxNElement:word read MaxFNElement write SetMaxNElement;
    {максимально возможная разница между узлами}
    property MaxWidth:word read MaxFWidth write SetMemory;
    {число элементов нагруженных распределенной нагрузкой}
    property NLoadElements:word read FNLoadElements write SetNLoadElements;
    {число нагруженных точек}
    property NLoadPoints:word read FNLoadPoints write SetNLoadPoints;
    {число узлов лежащих на пути}
    property NPPath:word read FNPPath write SetNPPath;
    {число Гауссовых точек в элементе}
    property NGaussPoints:byte read FNGaussPoints write SetNGaussPoints;
{---------------------------------------------------------}
    constructor Create(NPnt,NFr,NG:word);
    destructor Destroy;override;
    procedure GaussIntegSquare(var M:array of Extended;NEl:word;Func:TFindMatrix);
    procedure GaussIntegLine(var M:array of Extended;NEl:word;Func:TFindMatrix);
    procedure FormMatrixFind(var z:array of Extended;ksi,nu:Extended);{Матрицу формы}
    procedure Jacobian(var z:array of Extended; ksi,nu:Extended;x,y:array of Extended);
    {находит Якобиан}
    procedure DifFormGlob(var z:array of Extended;x,y:array of Extended;ksi,nu:Extended);{нахождение матрицы производной формы}
    procedure DeformationFind(var M:array of Extended;NEl:word;ksi,nu:Extended);{Нахождение вектора деформации}
    procedure RigidMatrixFind(var M:array of Extended;NEl:word;NGP:byte);{Нахождение матрицы жесткости в точке}
    procedure ApplyOnPoint(NPnt:word;NAxs,Sort:byte;Value:Extended);
    procedure CalcGlobal_B;
    {------------------ Virtual procedure ------------------}
    procedure FuncForm(var z:array of Extended;ksi,nu:Extended);virtual;{Функции формы}
    procedure DifFormLoc(var z:array of Extended;ksi,nu:Extended);virtual;{функция формы}
    procedure GukMatrixFind(var M:array of Extended;NEl:word;NGp:byte);virtual;
    procedure GradientMatrixFind(var M:array of Extended;NEl:word;ksi,nu:Extended);virtual;{Нахождение матрицы градиентов}
    procedure RigidMatrixAndForceFind(var M:array of Extended;NEl:word;NGp:byte);virtual;{Нахождение матрицы жесткости и вектора нагрузки}
    {----------------- Basic procedure ---------------------}
    procedure ElementFind(var M,V:array of Extended;NEl:word);dynamic;{Интегрирование по площади матрицы и вектора нагрузки для элемента}
    procedure FormatGlobalMatrix;dynamic;{Формирование глобальной матрицы жесткости}
    procedure Solve;dynamic;
    procedure OutDate(z:array of Extended;Stroka,NameFileOutDate:string);{Вывод данных}

  protected
    NPoint,NPoint2,l0,PointPlace,MaxNum,MinNum,DifNum:word;
    DxyDs,DksinuDs:array[0..1]of Extended;
    Delta:Extended;{определитель якобиана}
    Value0:array of string;
    ValueMatrix:array of Extended;
    Varabl,ksi,nu,dXPnt,dYPnt:Extended;
    FileValue:TextFile;{}
    Text:string;{}
    NumPoint:array of Longword;
    RigidEl:array of Extended;{Матрица жесткости элемента}
    ForceEl,VectorEl:array of Extended;{Вектор сил элемента}
  end;

implementation

constructor TMKE.Create(NPnt,NFr,NG:word);
begin
  N:=NPnt;
  NFree:=NFr;
  NGuk:=NG;
  SetLength(DifLoc,2*N);
  SetLength(Guk,NGuk*NGuk);
  SetLength(Matrix_B,NGuk*NFree*N);
  SetLength(Matrix_BT,NGuk*NFree*N);
  SetLength(Matrix_BE,NGuk*NFree*N);
  SetLength(ElemMatrix,N*NFree*(N*NFree+1));
  SetLength(NumPoint,NFree*N);
  SetLength(RigidEl,NFree*NFree*N*N);
  SetLength(VectorEl,NFree*N);
  SetLength(ForceEl,NFree*N);
  SetLength(FormVector,N);
  SetLength(xp,N);
  SetLength(yp,N);
  CellVolume:=10;
end;

destructor TMKE.Destroy;
begin
  MassivPoint:=Nil;
  ElementPoint:=Nil;
  ElementMaterial:=Nil;
  LoadElements:=Nil;
  LoadPoints:=Nil;
  Path:=Nil;
  if MatrixRigid<>Nil then MatrixRigid.Destroy;
  if Vector<>Nil then Vector.Destroy;
  if Force<> Nil then Force.Destroy;
  Global_B:=Nil;
  PntG:=Nil;
  WG:=Nil;
  xp:=Nil;
  yp:=Nil;
  DifLoc:=Nil;
  FormVector:=Nil;
  Guk:=Nil;
  Matrix_B:=Nil;
  Matrix_BT:=Nil;
  Matrix_BE:=Nil;
  ElemMatrix:=Nil;
  Value0:=Nil;
  ValueMatrix:=Nil;
  NumPoint:=Nil;
  RigidEl:=Nil;
  VectorEl:=Nil;
  ForceEl:=Nil;
end;

procedure TMKE.GaussIntegSquare(var M:array of Extended;NEl:word;
                                   Func:TFindMatrix); {need - xp[],yp[]}
var
  i,j:word;
begin
  SetLength(ValueMatrix,High(M)+1);
  for j:=0 to High(M) do M[j]:=0;
  for i:=0 to NGaussPoints-1 do
    begin
      Func(ValueMatrix,NEl,i);
      Jacobian(Jacob,PntG[2*i],PntG[2*i+1],xp,yp);
      Delta:=DefMatrix(Jacob,2);
      for j:=0 to High(M) do M[j]:=M[j]+WG[i]*Delta*ValueMatrix[j];
    end;
  ValueMatrix:=Nil;
end;

procedure TMKE.GaussIntegLine(var M:array of Extended;NEl:word;
                               Func:TFindMatrix);  {need - xp[],yp[]}
var
  i,j:word;
begin
  SetLength(ValueMatrix,High(M)+1);
  for j:=0 to High(M) do M[j]:=0;
  Func(ValueMatrix,NEl,0);
  Jacobian(Jacob,PntG[0],PntG[1],xp,yp);
  Transpon(Jacob,2,2,Inv_Jacob);
  for i:=0 to High(WG)-1 do
    begin
      DksinuDs[0]:=PntG[2*(i+1)]-PntG[i];
      DksinuDs[1]:=PntG[2*(i+1)+1]-PntG[i+1];
      MatrixSquareXVector(Inv_Jacob,DksinuDs,DxyDs);
      if i=0 then Delta:=WG[i]*LengthVector(DxyDs)
      else Delta:=0.5*WG[i]*LengthVector(DxyDs);
      for j:=0 to High(M) do M[j]:=M[j]+Delta*ValueMatrix[j];
      Func(ValueMatrix,NEl,i+1);
      Jacobian(Jacob,PntG[2*(i+1)],PntG[2*(i+1)+1],xp,yp);
      Transpon(Jacob,2,2,Inv_Jacob);
      MatrixSquareXVector(Inv_Jacob,DksinuDs,DxyDs);
      if i=(High(WG)-1) then Delta:=WG[i+1]*LengthVector(DxyDs)
      else Delta:=0.5*WG[i+1]*LengthVector(DxyDs);
      for j:=0 to High(M) do M[j]:=M[j]+Delta*ValueMatrix[j];
    end;
  ValueMatrix:=Nil;
end;

procedure TMKE.FormMatrixFind(var z:array of Extended;ksi,nu:Extended);
var
  i,j:word;
begin
  FuncForm(FormVector,ksi,nu);
  for i:=0 to NFree*NFree*N-1 do z[i]:=0;
  for i:=0 to N-1 do
    for j:=0 to NFree-1 do
      z[j*N*NFree+i*NFree+j]:=FormVector[i];
end;

procedure TMKE.Jacobian(var z:array of Extended;ksi,nu:Extended;x,y:array of Extended);
var
  i:integer;
begin
  DifFormLoc(DifLoc,ksi,nu);
  z[0]:=0;
  z[1]:=0;
  z[2]:=0;
  z[3]:=0;
  for i:=0 to N-1 do
    begin
      z[0]:=z[0]+Element(0,i,DifLoc,2,N)*x[i];
      z[1]:=z[1]+Element(0,i,DifLoc,2,N)*y[i];
      z[2]:=z[2]+Element(1,i,DifLoc,2,N)*x[i];
      z[3]:=z[3]+Element(1,i,DifLoc,2,N)*y[i];
    end;
end;

procedure TMKE.DifFormGlob(var z:array of Extended;x,y:array of Extended;
          ksi,nu:Extended);
begin
  Jacobian(Jacob,ksi,nu,x,y);
  OpositeMatrix(Jacob,Inv_Jacob,2);
  MatrixXMatrix(Inv_Jacob,DifLoc,2,2,N,z);
end;

procedure TMKE.RigidMatrixFind(var M:array of Extended;NEl:word;NGP:byte);
var
  i,j:byte;
begin
  GukMatrixFind(Guk,NEl,NGP);
  if Global_B=Nil then GradientMatrixFind(Matrix_B,NEl,PntG[2*NGP],PntG[2*NGP+1])
  else
    for i:=0 to NGuk-1 do
      for j:=0 to N*NFree-1 do Matrix_B[i*N*NFree+j]:=
       Global_B[NEl*FNGaussPoints*NGuk*N*NFree+NGP*NGuk*N*NFree+i*N*NFree+j];
  Transpon(Matrix_B,NGuk,NFree*N,Matrix_BT);
  MatrixXMatrix(Matrix_BT,Guk,NFree*N,NGuk,NGuk,Matrix_BE);
  MatrixXMatrix(Matrix_BE,Matrix_B,NFree*N,NGuk,NFree*N,M);
end;

procedure TMKE.ElementFind(var M,V:array of Extended;NEl:word);
var
  i,j:byte;
begin
  for i:=0 to N-1 do
    begin
      xp[i]:=MassivPoint[2*ElementPoint[N*NEl+i]];
      yp[i]:=MassivPoint[2*ElementPoint[N*NEl+i]+1];
    end;
  GaussIntegSquare(ElemMatrix,NEl,RigidMatrixAndForceFind);
  for i:=0 to N*NFree-1 do
    for j:=0 to N*NFree-1 do
      M[i*N*NFree+j]:=ElemMatrix[i*(N*NFree+1)+j];
  for i:=0 to NFree*N-1 do V[i]:=ElemMatrix[i*(N*NFree+1)+N*NFree];
end;

procedure TMKE.DeformationFind(var M:array of Extended;NEl:word;ksi,nu:Extended); {need - VectorEl[]}
begin
  GradientMatrixFind(Matrix_B,NEl,ksi,nu);
  MatrixXMatrix(Matrix_B,VectorEl,NGuk,NFree*N,1,M);
end;

procedure TMKE.ApplyOnPoint(NPnt:word;NAxs,Sort:byte;Value:Extended);
var
  l:integer;
begin
  if Sort=0 then Force.Massiv[NFree*Npnt+NAxs]:=Force.Massiv[NFree*NPnt+NAxs]+Value

  else if Sort=1 then
    begin
      Force.Massiv[NFree*Npnt+NAxs]:=Value;
      if not(MatrixRigid.Typ=LDL) then
        begin
          if  (NFree*NPnt+NAxs)<(MatrixRigid.Size2-1) then l0:=0
          else l0:=NFree*NPnt+NAxs-MatrixRigid.Size2+1;
          for l:=l0 to NFree*NPnt+NAxs-1 do
            begin
              Varabl:=MatrixRigid.Element[l,NFree*NPnt+NAxs];
              MatrixRigid.Element[l,NFree*NPnt+NAxs]:=0;
              Force.Massiv[l]:=Force.Massiv[l]-Varabl*Value;
            end;
          if (NFree*NPnt+NAxs)>(MatrixRigid.Size1-MatrixRigid.Size2)
             then l0:=MatrixRigid.Size1-1
          else l0:=NFree*NPnt+NAxs+MatrixRigid.Size2-1;
          for l:=NFree*NPnt+NAxs+1 to l0 do
            begin
              Varabl:=MatrixRigid.Element[NFree*NPnt+NAxs,l];
              MatrixRigid.Element[NFree*NPnt+NAxs,l]:=0;
              Force.Massiv[l]:=Force.Massiv[l]-Varabl*Value;
            end;
          MatrixRigid.Element[NFree*NPnt+NAxs,NFree*NPnt+NAxs]:=1;
        end;
    end

  else if (Sort=2) and (not(MatrixRigid.Typ=LDL)) then
    begin
      Varabl:=MatrixRigid.Element[NFree*Npnt+NAxs,NFree*NPnt+NAxs];
      MatrixRigid.Element[NFree*Npnt+NAxs,NFree*Npnt+NAxs]:=Varabl+Value;
    end;
end;

procedure TMKE.CalcGlobal_B;
var
  i:word;
  j,i0,j0:byte;
begin
  if Global_B=Nil then SetLength(Global_B,MaxFNElement*FNGaussPoints*NGuk*N*NFree);
  for i:=0 to NElement-1 do
    begin
      for j:=0 to N-1 do
        begin
          xp[j]:=MassivPoint[2*ElementPoint[i*N+j]];
          yp[j]:=MassivPoint[2*ElementPoint[i*N+j]+1];
        end;
      for j:=0 to FNGaussPoints-1 do
        begin
          GradientMatrixFind(Matrix_B,i,PntG[2*j],PntG[2*j+1]);
          for i0:=0 to NGuk-1 do
            for j0:=0 to N*NFree-1 do
              Global_B[i*FNGaussPoints*NGuk*N*NFree+j*NGuk*N*NFree+i0*N*NFree+j0]:=Matrix_B[i0*N*NFree+j0];
        end;
    end;
end;

{---------------------------------------------------------}
{-------------------property procedure--------------------}
{---------------------------------------------------------}
procedure TMKE.ReadingDate(Value:string);

begin
{-----------Чтение данных и подготовка к расчету---------}
  FNameFileDate:=Value;
  ProcessName:='Чтение данных';
  ProcessPos:=0;
  AssignFile(FileValue,Value);
  Reset(FileValue);
  ReadLn(FileValue,Text);
  if Pos('MyFormat',Text)>0 then ReadMyFormatDate;
  if Pos('Ansys',Text)>0 then ReadAnsysDate;
  CloseFile(FileValue);
  Value0:=Nil;
  ProcessName:='';
end;

procedure TMKE.ReadMyFormatDate;
var
  i,j:word;
begin
  ReadLn(FileValue,Text);
  repeat
    if Pos('CellVolume=',Text)>0 then CellVolume:=IntFromText(Text);
    if Pos('NumberOfPoints=',Text)>0 then MaxNAll:=IntFromText(Text);{Число узлов всего}
    if Pos('NumberOfElements=',Text)>0 then MaxNElement:=IntFromText(Text);{Число элементов}
    if Pos('MaxNumberDiff=',Text)>0 then MaxWidth:=IntFromText(Text);{Максимальная разница между номерами узлов}
    Readln(FileValue,Text);
  until Pos('[Points:]',Text)<>0;
  ProcessMaxPos:=NAll+NElement-1;
  ProcessPos:=1;
{-----------Загрузка координат узлов-----------}
  SetLength(Value0,2);
  for j:=0 to NAll-1 do
    begin
      ProcessPos:=1+j;
      Readln(FileValue,Text);
      StrokaValue(Text,';',Value0);
      for i:=0 to 1 do
        begin
          Varabl:=StrToFloat(Value0[i]);
          MassivPoint[j*2+i]:=Varabl;
          if (j=0)and(i=0) then Xmin:=Varabl
          else if (Xmin>Varabl)and(i=0) then Xmin:=Varabl;
          if (j=0)and(i=0) then Xmax:=Varabl
          else if (Xmax<Varabl)and(i=0) then Xmax:=Varabl;
          if (j=0)and(i=1) then Ymin:=Varabl
          else if (Ymin>Varabl)and(i=1) then Ymin:=Varabl;
          if (j=0)and(i=1) then Ymax:=Varabl
          else if (Ymax<Varabl)and(i=1) then Ymax:=Varabl;
        end;
    end;
  Value0:=Nil;
{----------Загрузка элементов------------------}
  repeat Readln(FileValue,Text) until Pos('[Elements:]',Text)<>0;
  SetLength(Value0,N+1);
  for j:=0 to NElement-1 do
    begin
      ProcessPos:=NAll+j;
      Readln(FileValue,Text);
      StrokaValue(Text,';',Value0);
      for i:=0 to N-1 do ElementPoint[j*N+i]:=StrToInt(Value0[i]);
      ElementMaterial[j]:=StrToInt(Value0[N]);
    end;
end;

procedure TMKE.ReadAnsysDate;
var
  i,j:word;
begin
  ReadLn(FileValue,Text);
  repeat
    if Pos('CellVolume=',Text)>0 then CellVolume:=IntFromText(Text);
    if Pos('NumberOfPoints=',Text)>0 then MaxNAll:=IntFromText(Text);{Число узлов всего}
    if Pos('NumberOfElements=',Text)>0 then MaxNElement:=IntFromText(Text);{Число элементов}
    Readln(FileValue,Text);
  until Pos('NODE',Text)=4;
  ProcessMaxPos:=NAll+NElement-1;
  ProcessPos:=1;
{-----------Загрузка координат узлов-----------}
  SetLength(Value0,2);
  j:=0;
  repeat
    Readln(FileValue,Text);
    if Pos(IntToStr(j+1),Text)<>0 then
      begin
        ProcessPos:=1+j;
        PointPlace:=Pos('.',Text);
        Text[PointPlace]:=',';
        PointPlace:=Pos('.',Text);
        Text[PointPlace]:=',';
        Value0[0]:=Copy(Text,12,15);
        Value0[1]:=Copy(Text,32,15);
        for i:=0 to 1 do
          begin
            Varabl:=StrToFloat(Value0[i]);
            MassivPoint[j*2+i]:=Varabl;
            if (j=0)and(i=0) then Xmin:=Varabl
            else if (Xmin>Varabl)and(i=0) then Xmin:=Varabl;
            if (j=0)and(i=0) then Xmax:=Varabl
            else if (Xmax<Varabl)and(i=0) then Xmax:=Varabl;
            if (j=0)and(i=1) then Ymin:=Varabl
            else if (Ymin>Varabl)and(i=1) then Ymin:=Varabl;
            if (j=0)and(i=1) then Ymax:=Varabl
            else if (Ymax<Varabl)and(i=1) then Ymax:=Varabl;
          end;
        j:=j+1;
      end;
  until Pos('ELEM',Text)=5;
  Value0:=Nil;
{----------Загрузка элементов------------------}
  SetLength(Value0,N+1);
  j:=0;
  repeat
    Readln(FileValue,Text);
    if Pos(IntToStr(j+1),Text)<>0 then
      begin
        ProcessPos:=NAll+j;
        Value0[N]:=Copy(Text,9,4);
        ElementMaterial[j]:=StrToInt(Value0[N])-1;
        Value0[0]:=Copy(Text,30,6);
        ElementPoint[j*N]:=StrToInt(Value0[0])-1;
        MaxNum:=ElementPoint[j*N];
        MinNum:=MaxNum;
        for i:=1 to N-1 do
          begin
            Value0[i]:=Copy(Text,30+6*i,6);
            ElementPoint[j*N+i]:=StrToInt(Value0[i])-1;
            if MaxNum<ElementPoint[j*N+i] then MaxNum:=ElementPoint[j*N+i];
            if MinNum>ElementPoint[j*N+i] then MinNum:=ElementPoint[j*N+i];
          end;
        if (j=0)or(DifNum<(MaxNum-MinNum)) then DifNum:=MaxNum-MinNum;
        j:=j+1;
      end;
  until j=MaxNElement;
  MaxWidth:=DifNum;{Максимальная разница между номерами узлов}
end;

procedure TMKE.FormatConcenLoad(Value:string);
var
  j:word;
begin
  {---------------------- ConcenLoad ----------------------}
  {Виды закрепления узлов: 0 - нет закрепления; 1 - X ; 2 - Y ; 3 - XY ;
  4 - w ; 5 - Xw ; 6 - Yw ; 7 - XYw}
  FNameFileConcen:=Value;
  AssignFile(FileValue,Value);
  Reset(FileValue);
  repeat
    Readln(FileValue,Text);
  until Pos('NumberOfPoints=',Text)>0;
  NLoadPoints:=IntFromText(Text);{Число элементов}
  ProcessName:='Формирование вектора нагрузки-2';
  ProcessPos:=0;
  ProcessMaxPos:=NLoadPoints;
  j:=0;
  while not EOF(FileValue)do
    begin
      ProcessPos:=j;
      Readln(FileValue,Text);
      if Pos('LoadType=1',Text)>0 then
        begin
          Readln(FileValue,Text);
          while Pos('End',Text)=0 do
            begin
              LoadPoints[j].Point:=StrToInt(Copy(Text,0,8))-1;
              LoadPoints[j].LType:=1;
              j:=j+1;
              Readln(FileValue,Text);
            end;
        end;
      if Pos('LoadType=2',Text)>0 then
        begin
          Readln(FileValue,Text);
          while Pos('End',Text)=0 do
            begin
              LoadPoints[j].Point:=StrToInt(Copy(Text,0,8))-1;
              LoadPoints[j].LType:=2;
              j:=j+1;
              Readln(FileValue,Text);
            end;
        end;
    end;
  CloseFile(FileValue);
  ProcessName:='';
end;

procedure TMKE.FormatDistribLoad(Value:string);
var
  i,j:word;
begin
  {---------------- Distribution Load -------------------}
  {Виды закрепления узлов: 0 - равномерное нормальное давление;1 - равномерное касательное давление}
  FNameFileDistrib:=Value;
  AssignFile(FileValue,Value);
  Reset(FileValue);
  ProcessName:='Формирование вектора нагрузки-1';
  ProcessPos:=0;
  Readln(FileValue,Text);
  repeat
    if Pos('NumberOfElements=',Text)>0 then NLoadElements:=IntFromText(Text);{Число элементов}
    Readln(FileValue,Text);
  until Pos('ELEM',Text)=5;
  ProcessMaxPos:=NLoadElements*2;
{----------Загрузка элементов------------------}
  j:=0;
  repeat
    Readln(FileValue,Text);
    if (Length(Text)>2)and(Pos('ELEM',Text)<>5) then
      begin
        ProcessPos:=j;
        LoadElements[(N+1)*j]:=StrToInt(Copy(Text,0,8))-1;
        j:=j+1;
      end;
  until j=NLoadElements;
{-----------Загрузка точек---------------------}
  j:=0;
  repeat
    Readln(FileValue,Text);
    if (Length(Text)>2)and(Pos('NODE',Text)<>4)and(Pos('LIST',Text)<>2)
      and(Pos('SORT',Text)<>2) then
      begin
        ProcessPos:=NLoadElements+j;
        for i:=0 to NLoadElements-1 do
          for j:=0 to N-1 do
            if (StrToInt(Copy(Text,0,8))-1)=ElementPoint[LoadElements[i*(N+1)]*N+j] then
              LoadElements[i*(N+1)+j+1]:=1;
      end;
  until EOF(FileValue);
end;

procedure TMKE.FormatPath(Value:string);
var
  i:word;
begin
  FNameFilePath:=Value;
  AssignFile(FileValue,Value);
  Reset(FileValue);
  ProcessName:='Формирование пути';
  ProcessPos:=0;
  Readln(FileValue,Text);
  NPPath:=IntFromText(Text);{Число узлов}
  ProcessMaxPos:=NPPath;
  for i:=0 to NPPath-1 do
    begin
      ProcessPos:=i;
      Readln(FileValue,Text);
      if Length(Text)<65 then Path[i]:=StrToInt(Text)
      else Path[i]:=StrToInt(Copy(Text,0,8))-1;
    end;
end;

procedure TMKE.SetNAll(Value:word);
begin
  FNAll:=Value;
  MatrixRigid.Size1:=NFree*FNAll;
  Vector.Size1:=NFree*FNAll;
  Force.Size1:=NFree*FNAll;
end;

procedure TMKE.SetWidth(Value:word);
begin
  FWidth:=Value;
  MatrixRigid.Size2:=NFree*(FWidth+1);
end;

procedure TMKE.SetMaxNAll(Value:word);
begin
  MaxFNAll:=Value;
  FNAll:=Value;
  SetLength(MassivPoint,2*Value);
  Vector:=TMatrix.Create(Usual,NFree*Value,1,CellVolume);
  Force:=TMatrix.Create(Usual,NFree*Value,1,CellVolume);
end;

procedure TMKE.SetMaxNElement(Value:word);
begin
  MaxFNElement:=Value;
  NElement:=Value;
  SetLength(ElementPoint,NElement*N);
  SetLength(ElementMaterial,NElement);
end;

procedure TMKE.SetNLoadElements(Value:word);
begin
  FNLoadElements:=Value;
  SetLength(LoadElements,(N+1)*Value);
end;

procedure TMKE.SetNLoadPoints(Value:word);
begin
  FNLoadPoints:=Value;
  SetLength(LoadPoints,Value);
end;

procedure TMKE.SetNPPath(Value:word);
begin
  Path:=Nil;
  FNPPath:=Value;
  SetLength(Path,Value);
end;

procedure TMKE.SetMemory(Value:word);
begin
  MaxFWidth:=Value;
  FWidth:=Value;
  MatrixRigid:=TMatrix.Create(Sym,NFree*MaxFNAll,NFree*(MaxFWidth+1),CellVolume);
end;

procedure TMKE.SetNGaussPoints(Value:byte);
begin
  FNGaussPoints:=Value;
  SetLength(PntG,2*Value);
  SetLength(WG,Value);
end;

{---------------------------------------------------------}
{-----------------------virtual----------------------------}
{---------------------------------------------------------}
procedure TMKE.FuncForm(var z:array of Extended;ksi,nu:Extended);
begin
{Функции формы}
end;

procedure TMKE.DifFormLoc(var z:array of Extended;ksi,nu:Extended);
begin
{Функция формы}
end;

procedure TMKE.GukMatrixFind(var M:array of Extended;NEl:word;NGp:byte);
begin
{Матрица констант}
end;

procedure TMKE.GradientMatrixFind(var M:array of Extended;NEl:word;ksi,nu:Extended);
begin
{Нахождение матрицы градиентов}
end;

procedure TMKE.RigidMatrixAndForceFind(var M:array of Extended;NEl:word;NGP:byte);
begin
{Нахождение матрицы жесткости и вектора нагрузки}
end;

{--------------------------------------------------------}
{---------------- основные процедуры --------------------}
{--------------------------------------------------------}
procedure TMKE.FormatGlobalMatrix;
var
  i,j,l:word;
begin
  ProcessName:='Формирование глобальной матрицы жесткости';
  ProcessPos:=0;
  ProcessMaxPos:=NFree*NAll+NElement-1;
  for i:=0 to NFree*NAll-1 do
    begin
      ProcessPos:=i;
      Force.Massiv[i]:=0;
    end;
  MatrixRigid.Zero;
  MatrixRigid.Typ:=Sym;
  for j:=0 to NElement-1 do
    begin
      ProcessPos:=NFree*NAll+j;
      {-------------------------------------------------------}
      {---------Нахождение матрицы жесткости элемента---------}
      ElementFind(RigidEl,ForceEl,j);
      {----------------------------------------------------------}
      {---------заполнение глобальной матрицы жесткости-----------------------}
      for i:=0 to N-1 do
        for l:=0 to NFree-1 do
          begin
            NPoint:=ElementPoint[j*N+i];
            NumPoint[i*NFree+l]:=NFree*NPoint+l;
            Force.Massiv[NFree*NPoint+l]:=Force.Massiv[NFree*NPoint+l]+ForceEl[NFree*i+l]
          end;
      MatrixRigid.AddSymBlock(RigidEl,NumPoint,NFree*N);
    end;
  ProcessName:='';
end;

procedure TMKE.Solve;
begin
  if MatrixRigid.Typ<>LDL then MatrixRigid.TransformToLDL;
  SolveLDL(MatrixRigid,Vector,Force);
end;

procedure TMKE.OutDate(z:array of Extended;Stroka,NameFileOutDate:string);
var
  i,j:word;
begin
  ProcessName:='Запись в файл';
  ProcessPos:=0;
  ProcessMaxPos:=2*N*NElement-N;
  Zmin:=z[0];
  Zmax:=z[0];
  for i:=1 to N*NElement-1 do
    begin
      ProcessPos:=i;
      if Zmax<z[i] then Zmax:=z[i];
      if Zmin>z[i] then Zmin:=z[i];
    end;
  AssignFile(FileValue,NameFileOutDate);
  Rewrite(FileValue);
  WriteLn(FileValue,Stroka);
  WriteLn(FileValue,'Число элементов:'+IntToStr(NElement));
  WriteLn(FileValue,'Число узлов в элементе:'+IntToStr(N));
  WriteLn(FileValue,'Минимальная координата по оси Х:'+FloatToStr(Xmin));
  WriteLn(FileValue,'Максимальная координата по оси Х:'+FloatToStr(Xmax));
  WriteLn(FileValue,'Минимальноя координата по оси У:'+FloatToStr(Ymin));
  WriteLn(FileValue,'Максимальная координата по оси У:'+FloatToStr(Ymax));
  WriteLn(FileValue,'Минимальное значение:'+FloatToStr(Zmin));
  WriteLn(FileValue,'Максимальное значение:'+FloatToStr(Zmax));
  for i:=0 to NElement-1 do
    begin
      ProcessPos:=N*NElement+N*i;
      for j:=0 to N-1 do
        begin
          NPoint:=ElementPoint[N*i+j];
          Write(FileValue,FloatToStr(MassivPoint[2*NPoint+0])+'|'+
               FloatToStr(MassivPoint[2*NPoint+1])+'|'+FloatToStr(z[i*N+j]));
          if j=(N-1) then WriteLn(FileValue,'|')
          else Write(FileValue,'|');
        end;
    end;
  CloseFile(FileValue);
end;

end.

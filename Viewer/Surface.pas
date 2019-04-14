unit Surface;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  Menus, ExtCtrls, ExtDlgs, ComCtrls, DdeMan;

type
  TSurfaceForm = class(TForm)
    FMainMenu: TMainMenu;
    FMenuFile: TMenuItem;
    FMenuProp: TMenuItem;
    FSaveFile: TMenuItem;
    FExit: TMenuItem;
    DrawingInMenu: TMenuItem;
    MenuColorSetka: TMenuItem;
    SetkaWhite: TMenuItem;
    SetkaBlack: TMenuItem;
    FParametr: TMenuItem;
    Open: TMenuItem;
    OpenForm: TOpenDialog;
    VisiblePictureMenu: TMenuItem;
    VisiblePalitraMenu: TMenuItem;
    VisibleSetkaMenu: TMenuItem;
    SavePicture: TSavePictureDialog;
    StatusBar3: TStatusBar;
    ScrollBox: TScrollBox;
    FImage: TImage;
    DdeServer: TDdeServerConv;

    {Аппроксимирующая функция}
    procedure Approx(var z0:array of Extended; ksi,nu:Extended);
     {Процедура прорисовки элемента}
    procedure DrawingElement(x,y,z:array of Extended);
    {Процедура прорисовки палитры}
    procedure DrawingPalitra(zmax,zmin:Extended);
    {Процедура прорисовки по данным находящимся в файле NameFile}
    procedure Drawing(FileName:string);

    {Событие - нажатие на "Выход"}
    procedure ExitClick(Sender:TObject);
    {Событик - в меню "Рисовать" нажатие на "Рисунок"}
    procedure VisiblePictureMenuClick(Sender:TObject);
    {Событие - в меню "Рисовать" нажатие на "Палитра"}
    procedure VisiblePalitraMenuClick(Sender:TObject);
    {Событие - в меню "Рисовать" нажатие на "Сетка" }
    procedure VisibleSetkaMenuClick(Sender:TObject);
    {Событие - в меню "Цвет сетки" нажатие на "Белый"}
    procedure SetkaWhiteClick(Sender: TObject);
    {Событие - в меню "Цвет сетки" нажатие на "Черный"}
    procedure SetkaBlackClick(Sender: TObject);
    {Событие - нажатие "Открыть"}
    procedure OpenClick(Sender: TObject);
    {Событие - нажатие "Сохранить"}
    procedure SaveClick(Sender: TObject);
    {Событие - нажатие на "Параметры"}
    procedure ParametrClick(Sender:TObject);
    procedure DdeServerExecuteMacro(Sender: TObject; Msg: TStrings);

  private
    { Private declarations }

    FVisiblePicture:Boolean;{Видимость рисунка}
    FVisiblePalitra:Boolean;{Видимость палитры}
    FVisibleSetka:Boolean;{Видимость сетки элементов}
    FColorSetka:TColor;{Цвет сетки элементов}

    {Установка свойства видимости рисунка}
    procedure SetVisiblePicture(Value:Boolean);
    {Установка свойства видимости палитры}
    procedure SetVisiblePalitra(Value:Boolean);
    {Установка свойства видимости сетки элементов}
    procedure SetVisibleSetka(Value:Boolean);
    {Установка свойства цвета сетки элементов}
    procedure SetColorSetka(Value:TColor);

  public
    { Public declarations }
    FFileNameDate:string;{Имя файла данных}
    FXPix:integer;{Ширина рисунка в пикселах}
    FYPix:integer;{Высота рисунка в пикселах}
    FAutoLimit:Boolean;{Авто определение пределов палитры}
    FZmax:Extended;{Максимальное значение палитры}
    FZmin:Extended;{Минимальное значение палитры}
    BrushKff:Extended;{Коэффициент заполнения элемента}
    MaxNumberColor:word;{Число цветов в палитре}
    ListColor: array[0..100] of TColor;{Палитра}

  protected
    { Protected declarations }
    FileDate:TextFile;
    Flag:Boolean;
    Area:TRect;
    Stroka,TextValue:string;
    NElement,n,HFont,HPalitra,Nnumber,Npoint,Naxis,
    Nksi,pix_x,pix_y,NColorPoint:integer;
    Scale,dS,dSmax,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,ValueZ:Extended;
    x,y,z:array of Extended;
    FuncForm:array of Extended;

  published
    { Published declarations }

    {Видимость рисунка}
    property VisiblePicture:Boolean read FVisiblePicture
             write SetVisiblePicture default True;
    {Видимость палитры}
    property VisiblePalitra:Boolean read FVisiblePalitra
             write SetVisiblePalitra default True;
    {Видимость сетки элементов}
    property VisibleSetka:Boolean read FVisibleSetka
             write SetVisibleSetka default False;
    {Цвет сетки}
    property ColorSetka:TColor read FColorSetka
             write SetColorSetka default clBlack;

  end;

const
  LengthExtended=9;
  Pix=20;

var
  SurfaceForm: TSurfaceForm;

implementation

uses Matrix_Rigid, UnitParametr;
{$R *.DFM}

procedure TSurfaceForm.Approx(var z0:array of Extended;ksi,nu:Extended);
begin
  if n=3 then
    begin
      z0[0]:=1-ksi-nu;
      z0[1]:=ksi;
      z0[2]:=nu;
    end;
  if n=4 then
    begin
      z0[0]:=(1-ksi)*(1-nu)/4;
      z0[1]:=(1+ksi)*(1-nu)/4;
      z0[2]:=(1+ksi)*(1+nu)/4;
      z0[3]:=(1-ksi)*(1+nu)/4;
    end;
end;

{Процедура прорисовки элемента}
procedure TSurfaceForm.DrawingElement(x,y,z:array of Extended);
var
  i,j,i0,j0,jk:integer;
begin
  if VisiblePicture then
    begin
      dSmax:=Sqrt((x[0]-x[N-1])*(x[0]-x[N-1])+(y[0]-y[N-1])*(y[0]-y[N-1]));
      for i:=0 to N-2 do
        begin
          dS:=Sqrt((x[i]-x[i+1])*(x[i]-x[i+1])+(y[i]-y[i+1])*(y[i]-y[i+1]));
          if dSmax<dS then dSmax:=dS;
        end;
      if n=3 then Nksi:=Round(BrushKff*dSmax*Scale)
      else Nksi:=Round(BrushKff*dSmax*Scale/2);
      if n=3 then i0:=0
      else i0:=-NKsi;
      for i:=i0 to Nksi do
        begin
          if n=3 then jk:=Nksi-i
          else jk:=Nksi;
          if n=3 then j0:=0
          else j0:=-Nksi;
          for j:=j0 to jk do
            begin
              Approx(FuncForm,i/Nksi,j/Nksi);
              ValueZ:=ScalarVector(FuncForm,z);
              pix_x:=Pix+Round((ScalarVector(FuncForm,x)-Xmin)*Scale);
              pix_y:=FYPix-Pix-Round((ScalarVector(FuncForm,y)-Ymin)*Scale);
              if FAutoLimit then NColorPoint:=
               Round((ValueZ-Zmin)*(MaxNumberColor-1)/(Zmax-Zmin))
              else
                begin
                  if ValueZ<FZmin then NColorPoint:=0
                  else if ValueZ>FZmax then NColorPoint:=MaxNumberColor-1
                  else NColorPoint:=Round((ValueZ-FZmin)*(MaxNumberColor-1)/(FZmax-FZmin));
                end;
              FImage.Picture.Bitmap.Canvas.Pixels[pix_x,pix_y]:=ListColor[NColorPoint];
            end;
        end;
    end;
  if VisibleSetka=True then
    begin
      FImage.Picture.Bitmap.Canvas.Pen.Color:=ColorSetka;
      pix_x:=Pix+Round((x[n-1]-Xmin)*Scale);
      pix_y:=FYPix-Pix-Round((y[n-1]-Ymin)*Scale);
      FImage.Picture.Bitmap.Canvas.MoveTo(pix_x,pix_y);
      for i:=0 to n-1 do
        begin
          pix_x:=Pix+Round((x[i]-Xmin)*Scale);
          pix_y:=FYPix-Pix-Round((y[i]-Ymin)*Scale);
          FImage.Picture.Bitmap.Canvas.LineTo(pix_x,pix_y);
        end;
    end;
end;

{Процедура прорисовки палитры}
procedure TSurfaceForm.DrawingPalitra(zmax,zmin:Extended);
var
  j:integer;
begin
  HPalitra:=FYPix-4*Pix;
  Nnumber:=HPalitra div (5*8);

  FImage.Picture.Bitmap.Canvas.Pen.Style:=psSolid;
  FImage.Picture.Bitmap.Canvas.Pen.Width:=1;
  FImage.Picture.Bitmap.Canvas.Brush.Style:=bsSolid;
  for j:=1 to MaxNumberColor do
    begin
      FImage.Picture.Bitmap.Canvas.Brush.Color:=ListColor[j-1];
      FImage.Picture.Bitmap.Canvas.Pen.Color:=ListColor[j-1];
      FImage.Picture.Bitmap.Canvas.Rectangle(FXPix-Pix-LengthExtended*8,FYPix-Pix-round(j*HPalitra/MaxNumberColor),
      FXPix-LengthExtended*8,FYPix-Pix-round((j-1)*HPalitra/MaxNumberColor));
    end;
  FImage.Picture.Bitmap.Canvas.Font.Height:=8;
  FImage.Picture.Bitmap.Canvas.Brush.Color:=clWhite;
  FImage.Picture.Bitmap.Canvas.Font.Color:=clBlack;
  for j:=0 to Nnumber do FImage.Picture.Bitmap.Canvas.TextOut(FXPix-LengthExtended*8,
    FYPix-(round(j*HPalitra/Nnumber)+10)-Pix,
    '.  '+FloatToStrF((zmax-zmin)*j/Nnumber+zmin,ffGeneral,4,3));
end;


{Процедура прорисовки}
procedure TSurfaceForm.Drawing(FileName:string);
var
  i,j:integer;
begin
  StatusBar3.Panels[2].Text:='Идет прорисовка...';
  Area.Left:=0;
  Area.Top:=0;
  Area.Right:=FXPix;
  Area.Bottom:=FYPix;
  FImage.Width:=FXPix;
  FImage.Height:=FYPix;
  FImage.Picture.Bitmap.Height:=FYPix;
  FImage.Picture.Bitmap.Width:=FXPix;
  FImage.Picture.Bitmap.Canvas.Pen.Style:=psSolid;
  FImage.Picture.Bitmap.Canvas.Pen.Color:=clWhite;
  FImage.Picture.Bitmap.Canvas.Brush.Color:=clWhite;
  FImage.Picture.Bitmap.Canvas.Brush.Style:=bsSolid;
  FImage.Picture.Bitmap.Canvas.Rectangle(0,0,FXPix,FYPix);

  AssignFile(FileDate,FileName);
  Reset(FileDate);
  ReadLn(FileDate,Stroka);
  HFont:=3*(FXPix-2*Pix) div 2*Length(Stroka);
  if HFont>Pix then HFont:=Pix;
  if HFont<8 then HFont:=8;
  FImage.Picture.Bitmap.Canvas.Font.Color:=clBlack;
  FImage.Picture.Bitmap.Canvas.Font.Height:=HFont;
  FImage.Picture.Bitmap.Canvas.TextOut(Pix,Pix,Stroka);

  ReadLn(FileDate,Stroka);
  TextValue:='';
  for i:=1 to Length(Stroka) do
    if (Stroka[i]<='9')and(Stroka[i]>='0') then TextValue:=TextValue+Stroka[i];
  NElement:=StrToInt(TextValue);

  ReadLn(FileDate,Stroka);
  TextValue:='';
  for i:=1 to Length(Stroka) do
    if (Stroka[i]<='9')and(Stroka[i]>='0') then TextValue:=TextValue+Stroka[i];
  n:=StrToInt(TextValue);
  SetLength(x,n);
  SetLength(y,n);
  SetLength(z,n);
  SetLength(FuncForm,n);

  for j:=0 to 5 do
    begin
      ReadLn(FileDate,Stroka);
      TextValue:='';
      Flag:=False;
      for i:=1 to Length(stroka) do
        begin
         if Flag and (Stroka<>' ') then TextValue:=TextValue+Stroka[i];
         if stroka[i]=':' then Flag:=True;
        end;
      case j of
        0: Xmin:=StrToFloat(TextValue);
        1: Xmax:=StrToFloat(TextValue);
        2: Ymin:=StrToFloat(TextValue);
        3: Ymax:=StrToFloat(TextValue);
        4: Zmin:=StrToFloat(TextValue);
        5: Zmax:=StrToFloat(TextValue);
      end;
    end;

  if VisiblePalitra=True then
    if FAutoLimit=True then DrawingPalitra(Zmax,Zmin)
    else if FAutoLimit=False then DrawingPalitra(FZmax,FZmin);

  Scale:=(FXPix-3*Pix-8*LengthExtended)/(Xmax-Xmin);
  if Scale>(FYPix-4*Pix)/(Ymax-Ymin) then Scale:=(FYPix-4*Pix)/(Ymax-Ymin);

  StatusBar3.Canvas.Pen.Color:=clNavy;
  StatusBar3.Canvas.Brush.Color:=clNavy;
  for j:=1 to NElement do
    begin
      StatusBar3.Canvas.Rectangle(156,6,157+Round(190*j/NElement),16);
      Application.ProcessMessages;
      ReadLn(FileDate,stroka);
      TextValue:='';
      NPoint:=0;
      Naxis:=0;
      for i:=1 to Length(stroka) do
        if stroka[i]='|' then
         begin
           case Naxis of
             0: x[NPoint]:=StrToFloat(TextValue);
             1: y[NPoint]:=StrToFloat(TextValue);
             2: z[NPoint]:=StrToFloat(TextValue);
           end;
           TextValue:='';
           if Naxis=2 then NPoint:=NPoint+1;
           if Naxis<2 then Naxis:=Naxis+1
           else Naxis:=0;
         end
        else if stroka[i]<>' ' then TextValue:=TextValue+stroka[i];
      DrawingElement(x,y,z);
    end;
  x:=Nil;
  y:=Nil;
  z:=Nil;
  FuncForm:=Nil;
  CloseFile(FileDate);
  StatusBar3.Panels[2].Text:='';
  StatusBar3.Canvas.Pen.Color:=clBtnFace;
  StatusBar3.Canvas.Brush.Color:=clBtnFace;
  StatusBar3.Canvas.Rectangle(156,6,347,16);
end;

{-----------------------------------------------}
{------ ПРОЦЕДУРЫ УСТАНОВКИ СВОЙСТВ ------------}
{-----------------------------------------------}
{Установка свойства видимости рисунка}
procedure TSurfaceForm.SetVisiblePicture(Value:Boolean);
begin
  if FVisiblePicture<>Value then
  begin
    FVisiblePicture:=Value;
    if (Value=False)and(VisibleSetka=False)then VisibleSetka:=True;
    if Value=False then VisiblePalitra:=False;
    VisiblePictureMenu.Checked:=Value;
  end;
end;

{Установка свойства видимости палитры}
procedure TSurfaceForm.SetVisiblePalitra(Value:Boolean);
begin
  if FVisiblePalitra<>Value then
  begin
    FVisiblePalitra:=Value;
    if (Value=True)and(VisiblePicture=False)then VisiblePicture:=True;
    VisiblePalitraMenu.Checked:=Value;
  end;
end;

{Установка свойства видимости сетки элементов}
procedure TSurfaceForm.SetVisibleSetka(Value:Boolean);
begin
  if FVisibleSetka<>Value then
  begin
    FVisibleSetka:=Value;
    if (Value=False)and(VisiblePicture=False)then VisiblePicture:=True;
    VisibleSetkaMenu.Checked:=Value;
  end;
end;

{Установка цвета сетки элементов}
procedure TSurfaceForm.SetColorSetka(Value:TColor);
begin
  if FColorSetka<>Value then
  begin
    FColorSetka:=Value;
    SetkaWhite.Checked:=False;
    SetkaBlack.Checked:=False;
    if Value=clBlack then SetkaBlack.Checked:=True
    else SetkaWhite.Checked:=True;
  end;
end;


{-------------------------------------------------}
{---- ПРОЦЕДУРЫ ПО ОБРАБОТКЕ СОБЫТИЙ -------------}
{-------------------------------------------------}
{нажатие "Выход"}
procedure TSurfaceForm.ExitClick(Sender:TObject);
begin
  Close;
end;

{в меню "Рисовать" нажатие на "Рисунок"}
procedure TSurfaceForm.VisiblePictureMenuClick(Sender:TObject);
begin
  VisiblePicture:=not(VisiblePicture);
  if FFileNameDate<>'' then Drawing(FFileNameDate);
end;

{в меню "Рисовать" нажатие на "Палитра"}
procedure TSurfaceForm.VisiblePalitraMenuClick(Sender:TObject);
begin
  VisiblePalitra:=not(VisiblePalitra);
  if FFileNameDate<>'' then Drawing(FFileNameDate);
end;

{в меню "Рисовать" нажатие на "Сетка" }
procedure TSurfaceForm.VisibleSetkaMenuClick(Sender:TObject);
begin
  VisibleSetka:=not(VisibleSetka);
  if FFileNameDate<>'' then Drawing(FFileNameDate);
end;

{в меню "Цвет сетки" нажатие на "Белый"}
procedure TSurfaceForm.SetkaWhiteClick(Sender: TObject);
begin
  ColorSetka:=clWhite;
end;

{в меню "Цвет сетки" нажатие на "Черный"}
procedure TSurfaceForm.SetkaBlackClick(Sender: TObject);
begin
  ColorSetka:=clBlack;
end;

{Нажатие на "Открыть"}
procedure TSurfaceForm.OpenClick(Sender: TObject);
begin
  if OpenForm.Execute then
   begin
     FFileNameDate:=OpenForm.FileName;
     OpenForm.InitialDir:=OpenForm.FileName;
     FSaveFile.Enabled:=True;
     Drawing(OpenForm.FileName);
   end;
end;

{нажатие "Сохранить"}
procedure TSurfaceForm.SaveClick(Sender: TObject);
begin
  if SavePicture.Execute then
    FImage.Picture.SaveToFile(SavePicture.FileName+'.bmp');
end;

{Событие - нажатие на "Параметры"}
procedure TSurfaceForm.ParametrClick(Sender:TObject);
begin
  FormParametr.Visible:=True;
end;

procedure TSurfaceForm.DdeServerExecuteMacro(Sender: TObject;
  Msg: TStrings);
begin
  {-------}
end;

end.

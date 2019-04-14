unit UnitParametr;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls, ExtCtrls, ComCtrls;

type
  TFormParametr = class(TForm)
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    ButtonApply: TButton;
    ButtonExit: TButton;
    UpDown1: TUpDown;
    Bevel: TBevel;
    CheckBox: TCheckBox;
    UpDown2: TUpDown;
    Edit1: TEdit;
    Edit2: TEdit;
    Edit4: TEdit;
    Edit5: TEdit;
    Label4: TLabel;
    Label5: TLabel;
    ComboBox1: TComboBox;
    Edit3: TEdit;
    Label6: TLabel;
    constructor Create(AOwner:TComponent);override;
    procedure Change1(Sender:TObject;var AllowChange:Boolean);
    procedure Change2(Sender:TObject;var AllowChange:Boolean);
    procedure ClickCheckBox(Sender:TObject);
    procedure ClickExit(Sender:TObject);
    procedure ClickApply(Sender:TObject);

  private
    { Private declarations }
    FT:TextFile;
    Text,Dir:string;
    Indx,LText:word;

  public
    { Public declarations }

  published
    { Published declarations }

  end;

var
  FormParametr: TFormParametr;

implementation

{$R *.DFM}

uses Surface;

procedure TFormParametr.Change1(Sender:TObject;var AllowChange:Boolean);
begin
  if UpDown1.Position<10 then UpDown1.Increment:=1
  else if UpDown1.Position<100 then UpDown1.Increment:=5
  else UpDown1.Increment:=10;
  Edit1.Text:=IntToStr(UpDown1.Position);
end;

procedure TFormParametr.Change2(Sender:TObject;var AllowChange:Boolean);
begin
  if UpDown2.Position<10 then UpDown2.Increment:=1
  else if UpDown2.Position<100 then UpDown2.Increment:=5
  else UpDown2.Increment:=10;
  Edit2.Text:=IntToStr(UpDown2.Position);
end;

procedure TFormParametr.ClickCheckBox(Sender:TObject);
begin
  if CheckBox.Checked=False then
      begin
        Edit4.ReadOnly:=False;
        Edit4.Color:=clWindow;
        Edit5.ReadOnly:=False;
        Edit5.Color:=clWindow;
      end
    else
      begin
        Edit4.ReadOnly:=True;
        Edit4.Color:=clMenu;
        Edit5.ReadOnly:=True;
        Edit5.Color:=clMenu;
      end;
end;

procedure TFormParametr.ClickExit(Sender:TObject);
begin
  Visible:=False;
end;

constructor TFormParametr.Create(AOwner:TComponent);
begin
  inherited Create(AOwner);
  Dir:=GetCurrentDir;
  CheckBox.Checked:=True;
  UpDown1.Position:=700;
  Edit1.Text:='700';
  UpDown2.Position:=500;
  Edit2.Text:='500';
  Edit3.Text:='1,4';
  Edit4.Text:='0';
  Edit5.Text:='0';
  AssignFile(FT,'Colors.dat');
  Reset(FT);
  while not EOF(FT) do
    begin
      ReadLn(FT,Text);
      Indx:=Pos('=',Text);
      if Indx<>0 then ComboBox1.Items.Add(Copy(Text,0,Indx-1));
    end;
  CloseFile(FT);
  ComboBox1.Text:=ComboBox1.Items.Strings[0];
  ClickApply(ButtonApply);
  with SurfaceForm do
    begin
      FFileNameDate:='';
      VisiblePicture:=True;
      VisibleSetka:=False;
      VisiblePalitra:=True;
      ColorSetka:=clBlack;
    end;
end;

procedure TFormParametr.ClickApply(Sender:TObject);
var
  i:word;
begin
  SurfaceForm.FAutoLimit:=CheckBox.Checked;
  SurfaceForm.FXPix:=UpDown1.Position;
  SurfaceForm.FYPix:=UpDown2.Position;
  SurfaceForm.BrushKff:=StrToFloat(Edit3.Text);
  SurfaceForm.FZmax:=StrToFloat(Edit4.Text);
  SurfaceForm.FZmin:=StrToFloat(Edit5.Text);
  AssignFile(FT,Dir+'\Colors.dat');
  Reset(FT);
  repeat
    ReadLn(FT,Text);
  until Pos(ComboBox1.Text,Text)<>0;
  Indx:=Pos('=',Text);
  LText:=Length(Text);
  Text:=Copy(Text,Indx+1,LText-Indx);
  Indx:=StrToInt(Text);
  SurfaceForm.MaxNumberColor:=Indx;
  for i:=0 to Indx-1 do
    begin
      ReadLn(FT,Text);
      SurfaceForm.ListColor[i]:=StrToInt(Text);
    end;
  CloseFile(FT);
  Visible:=False;
  if FileExists(SurfaceForm.FFileNameDate) then SurfaceForm.Drawing(SurfaceForm.FFileNameDate);
end;

end.

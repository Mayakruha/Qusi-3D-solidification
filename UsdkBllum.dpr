program UsdkBllum;

uses
  Forms,
  MainForm in 'MainForm.pas' {FormCheck};

{$R *.RES}

begin
  Application.Initialize;
  Application.CreateForm(TFormCheck, FormCheck);
  Application.Run;
end.

program ProjectSurface;

uses
  Forms,
  Surface in 'Surface.pas' {SurfaceForm},
  UnitParametr in 'UnitParametr.pas' {FormParametr};

{$R *.RES}

begin
  Application.Initialize;
  Application.CreateForm(TSurfaceForm, SurfaceForm);
  Application.CreateForm(TFormParametr, FormParametr);
  Application.Run;
end.

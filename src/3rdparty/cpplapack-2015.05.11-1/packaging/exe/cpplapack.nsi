!define version "2010.03.27-3"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;; base setting ;;;;;;;;
!include MUI2.nsh
Name "CPPLapack"
OutFile "cpplapack-${version}-setup.exe"
InstallDir "c:\cpplapack"

;;;;;;;; welcome ;;;;;;;;
;!define MUI_WELCOMEFINISHPAGE_BITMAP "C:\Program Files\NSIS\Contrib\Graphics\Wizard\win.bmp"

!define MUI_WELCOMEFINISHPAGE_BITMAP "${NSISDIR}\Contrib\Graphics\Wizard\win.bmp"
!define MUI_WELCOMEPAGE_TITLE "Welcome to CPPLapack installer."
!define MUI_WELCOMEPAGE_TEXT "Press $\"next$\" to procced."
!insertmacro MUI_PAGE_WELCOME

;;!insertmacro MUI_PAGE_LICENSE ".\LICENSE.txt"
;;!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_DIRECTORY
;;!insertmacro MUI_PAGE_STARTMENU
!insertmacro MUI_PAGE_INSTFILES

;;;;;;;; finish ;;;;;;;;
!define MUI_FINISHPAGE_TITLE "CPPLapack was installed successfully."
!define MUI_FINISHPAGE_TEXT "Press $\"Finish$\"."
!define MUI_FINISHPAGE_LINK "Click here to see the online documentation."
!define MUI_FINISHPAGE_LINK_LOCATION "http://cpplapack.sourceforge.net/"
!define MUI_FINISHPAGE_SHOWREADME "$INSTDIR/README.txt"
!insertmacro MUI_PAGE_FINISH

;;;;;;;; language ;;;;;;;;
!insertmacro MUI_LANGUAGE "English"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
Section "Install"
  SetOutPath $INSTDIR
  File .\README.txt
  
  SetOutPath $INSTDIR\lib
  File /r .\win32\*
  
  SetOutPath $INSTDIR\lib64
  File /r .\win64\*

  SetOutPath $INSTDIR\vc_sample
!system "mkdir vc_sample"
!system "xcopy ..\..\test\vc_sample .\vc_sample /s /e /q /r /y"
  File /r .\vc_sample\*
!system "rmdir /s /q .\vc_sample"

  SetOutPath $INSTDIR\include
  File .\stdint.h
  
  SetOutPath $INSTDIR
!system "mkdir include"
!system "xcopy ..\..\include .\include /s /e /q /r /y"
  File /r .\include
!system "rmdir /s /q include"

;;  WriteUninstaller $INSTDIR\uninstall.exe
SectionEnd

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Section "Uninstall"
;  Delete $INSTDIR\uninstall.exe
;  RMDir /r $INSTDIR
;SectionEnd

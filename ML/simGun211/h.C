// Mainframe macro generated from application: /root/alice/sw/slc8_x86-64/ROOT/v6-28-02-7/bin/root.exe
// By ROOT version 6.28/02 on 2023-07-06 11:23:23

#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif
#ifndef ROOT_TGListBox
#include "TGListBox.h"
#endif
#ifndef ROOT_TGScrollBar
#include "TGScrollBar.h"
#endif
#ifndef ROOT_TGComboBox
#include "TGComboBox.h"
#endif
#ifndef ROOT_TGMenu
#include "TGMenu.h"
#endif
#ifndef ROOT_TGFileDialog
#include "TGFileDialog.h"
#endif
#ifndef ROOT_TGButtonGroup
#include "TGButtonGroup.h"
#endif
#ifndef ROOT_TGCanvas
#include "TGCanvas.h"
#endif
#ifndef ROOT_TGFSContainer
#include "TGFSContainer.h"
#endif
#ifndef ROOT_TGButton
#include "TGButton.h"
#endif
#ifndef ROOT_TRootContextMenu
#include "TRootContextMenu.h"
#endif
#ifndef ROOT_TGFSComboBox
#include "TGFSComboBox.h"
#endif
#ifndef ROOT_TGLabel
#include "TGLabel.h"
#endif
#ifndef ROOT_TGListView
#include "TGListView.h"
#endif
#ifndef ROOT_TGSplitter
#include "TGSplitter.h"
#endif
#ifndef ROOT_TGTextEntry
#include "TGTextEntry.h"
#endif
#ifndef ROOT_TRootCanvas
#include "TRootCanvas.h"
#endif
#ifndef ROOT_TGDockableFrame
#include "TGDockableFrame.h"
#endif
#ifndef ROOT_TG3DLine
#include "TG3DLine.h"
#endif
#ifndef ROOT_TGToolTip
#include "TGToolTip.h"
#endif

#include "Riostream.h"

void h()
{

   // main frame
   TGMainFrame *fRootCanvas284 = new TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame);

   // menu bar
   TGMenuBar *fMenuBar294 = new TGMenuBar(fRootCanvas284,800,20,kHorizontalFrame);

   // "File" menu
   TGPopupMenu *fPopupMenu286 = new TGPopupMenu(gClient->GetDefaultRoot(),110,130,kOwnBackground);
   fPopupMenu286->AddEntry("&New Canvas",0);
   fPopupMenu286->AddEntry("&Open...",1);
   fPopupMenu286->AddEntry("&Close Canvas",13);
   fPopupMenu286->AddSeparator();

   // cascaded menu "Save"
   TGPopupMenu *fPopupMenu285 = new TGPopupMenu(gClient->GetDefaultRoot(),118,150,kOwnBackground);
   fPopupMenu285->AddEntry("tcnvRane0.&ps",5);
   fPopupMenu285->AddEntry("tcnvRane0.&eps",6);
   fPopupMenu285->AddEntry("tcnvRane0.p&df",7);
   fPopupMenu285->AddEntry("tcnvRane0.&tex",11);
   fPopupMenu285->AddEntry("tcnvRane0.&gif",8);
   fPopupMenu285->AddEntry("tcnvRane0.&jpg",9);
   fPopupMenu285->AddEntry("tcnvRane0.&png",10);
   fPopupMenu285->AddEntry("tcnvRane0.&C",4);
   fPopupMenu285->AddEntry("tcnvRane0.&root",3);
   fPopupMenu286->AddPopup("&Save",fPopupMenu285);
   fPopupMenu286->AddEntry("Save &As...",2);
   fPopupMenu286->AddSeparator();
   fPopupMenu286->AddEntry("&Print...",12);
   fPopupMenu286->AddSeparator();
   fPopupMenu286->AddEntry("&Quit ROOT",14);
   fMenuBar294->AddPopup("&File",fPopupMenu286, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "Edit" menu
   TGPopupMenu *fPopupMenu288 = new TGPopupMenu(gClient->GetDefaultRoot(),69,130,kOwnBackground);
   fPopupMenu288->AddEntry("&Style...",15);
   fPopupMenu288->AddSeparator();
   fPopupMenu288->AddEntry("Cu&t",16);
   fPopupMenu288->DisableEntry(16);
   fPopupMenu288->AddEntry("&Copy",17);
   fPopupMenu288->DisableEntry(17);
   fPopupMenu288->AddEntry("&Paste",18);
   fPopupMenu288->DisableEntry(18);
   fPopupMenu288->AddSeparator();

   // cascaded menu "Clear"
   TGPopupMenu *fPopupMenu287 = new TGPopupMenu(gClient->GetDefaultRoot(),74,38,kOwnBackground);
   fPopupMenu287->AddEntry("&Pad",19);
   fPopupMenu287->AddEntry("&Canvas",20);
   fPopupMenu288->AddPopup("C&lear",fPopupMenu287);
   fPopupMenu288->AddSeparator();
   fPopupMenu288->AddEntry("&Undo",21);
   fPopupMenu288->DisableEntry(21);
   fPopupMenu288->AddEntry("&Redo",22);
   fPopupMenu288->DisableEntry(22);
   fMenuBar294->AddPopup("&Edit",fPopupMenu288, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "View" menu
   TGPopupMenu *fPopupMenu290 = new TGPopupMenu(gClient->GetDefaultRoot(),122,162,kOwnBackground);
   fPopupMenu290->AddEntry("&Editor",23);
   fPopupMenu290->AddEntry("&Toolbar",24);
   fPopupMenu290->AddEntry("Event &Statusbar",25);
   fPopupMenu290->AddEntry("T&oolTip Info",26);
   fPopupMenu290->AddSeparator();
   fPopupMenu290->AddEntry("&Colors",27);
   fPopupMenu290->AddEntry("&Fonts",28);
   fPopupMenu290->DisableEntry(28);
   fPopupMenu290->AddEntry("&Markers",29);
   fPopupMenu290->AddSeparator();
   fPopupMenu290->AddEntry("&Iconify",30);
   fPopupMenu290->AddSeparator();

   // cascaded menu "View With"
   TGPopupMenu *fPopupMenu289 = new TGPopupMenu(gClient->GetDefaultRoot(),76,38,kOwnBackground);
   fPopupMenu289->AddEntry("&X3D",31);
   fPopupMenu289->AddEntry("&OpenGL",32);
   fPopupMenu290->AddPopup("&View With",fPopupMenu289);
   fMenuBar294->AddPopup("&View",fPopupMenu290, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "Options" menu
   TGPopupMenu *fPopupMenu291 = new TGPopupMenu(gClient->GetDefaultRoot(),148,194,kOwnBackground);
   fPopupMenu291->AddEntry("&Auto Resize Canvas",33);
   fPopupMenu291->CheckEntry(33);
   fPopupMenu291->AddEntry("&Resize Canvas",34);
   fPopupMenu291->AddEntry("&Move Opaque",35);
   fPopupMenu291->CheckEntry(35);
   fPopupMenu291->AddEntry("Resize &Opaque",36);
   fPopupMenu291->CheckEntry(36);
   fPopupMenu291->AddSeparator();
   fPopupMenu291->AddEntry("&Interrupt",37);
   fPopupMenu291->AddEntry("R&efresh",38);
   fPopupMenu291->AddSeparator();
   fPopupMenu291->AddEntry("&Pad Auto Exec",39);
   fPopupMenu291->AddSeparator();
   fPopupMenu291->AddEntry("&Statistics",40);
   fPopupMenu291->CheckEntry(40);
   fPopupMenu291->AddEntry("Histogram &Title",41);
   fPopupMenu291->CheckEntry(41);
   fPopupMenu291->AddEntry("&Fit Parameters",42);
   fPopupMenu291->AddEntry("Can Edit &Histograms",43);
   fMenuBar294->AddPopup("&Options",fPopupMenu291, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "Tools" menu
   TGPopupMenu *fPopupMenu292 = new TGPopupMenu(gClient->GetDefaultRoot(),120,102,kOwnBackground);
   fPopupMenu292->AddEntry("&Inspect ROOT",44);
   fPopupMenu292->AddEntry("&Class Tree",45);
   fPopupMenu292->AddEntry("&Fit Panel",46);
   fPopupMenu292->AddEntry("&Start Browser",47);
   fPopupMenu292->AddEntry("&Gui Builder",48);
   fPopupMenu292->AddEntry("&Event Recorder",49);
   fMenuBar294->AddPopup("&Tools",fPopupMenu292, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "Help" menu
   TGPopupMenu *fPopupMenu293 = new TGPopupMenu(gClient->GetDefaultRoot(),120,142,kOwnBackground);
   fPopupMenu293->AddLabel("Basic Help On...");
   fPopupMenu293->DefaultEntry(-1);
   fPopupMenu293->AddSeparator();
   fPopupMenu293->AddEntry("&Canvas",51);
   fPopupMenu293->AddEntry("&Menus",52);
   fPopupMenu293->AddEntry("&Graphics Editor",53);
   fPopupMenu293->AddEntry("&Browser",54);
   fPopupMenu293->AddEntry("&Objects",55);
   fPopupMenu293->AddEntry("&PostScript",56);
   fPopupMenu293->AddSeparator();
   fPopupMenu293->AddEntry("&About ROOT...",50);
   fMenuBar294->AddPopup("&Help",fPopupMenu293, new TGLayoutHints(kLHintsRight | kLHintsTop));
   fRootCanvas284->AddFrame(fMenuBar294, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,0,0,1,1));
   TGHorizontal3DLine *fHorizontal3DLine302 = new TGHorizontal3DLine(fRootCanvas284,1600,2);
   fRootCanvas284->AddFrame(fHorizontal3DLine302, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   // dockable frame
   TGDockableFrame *fDockableFrame303 = new TGDockableFrame(fRootCanvas284);

   // next lines belong to the dockable frame widget
   fDockableFrame303->EnableUndock(kTRUE);
   fDockableFrame303->EnableHide(kFALSE);
   fDockableFrame303->SetWindowName("ToolBar: tcnvRane0");
   fDockableFrame303->DockContainer();

   fRootCanvas284->AddFrame(fDockableFrame303, new TGLayoutHints(kLHintsExpandX));
   TGHorizontal3DLine *fHorizontal3DLine308 = new TGHorizontal3DLine(fRootCanvas284,1600,2);
   fRootCanvas284->AddFrame(fHorizontal3DLine308, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   // composite frame
   TGCompositeFrame *fCompositeFrame309 = new TGCompositeFrame(fRootCanvas284,800,514,kHorizontalFrame);

   // composite frame
   TGCompositeFrame *fCompositeFrame310 = new TGCompositeFrame(fCompositeFrame309,175,778,kFixedWidth);
   fCompositeFrame310->SetLayoutManager(new TGVerticalLayout(fCompositeFrame310));

   fCompositeFrame309->AddFrame(fCompositeFrame310, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));

   // canvas widget
   TGCanvas *fCanvas311 = new TGCanvas(fCompositeFrame309,800,514,kSunkenFrame);

   // canvas viewport
   TGViewPort *fViewPort312 = fCanvas311->GetViewPort();

   // canvas container
   Int_t canvasID = gVirtualX->InitWindow((ULongptr_t)fCanvas311->GetId());
   Window_t winC = gVirtualX->GetWindowID(canvasID);
   TGCompositeFrame *fCompositeFrame321 = new TGCompositeFrame(gClient,winC,fViewPort312);
   fViewPort312->AddFrame(fCompositeFrame321);
   fCompositeFrame321->SetLayoutManager(new TGVerticalLayout(fCompositeFrame321));
   fCompositeFrame321->MapSubwindows();
   fCanvas311->SetContainer(fCompositeFrame321);
   fCanvas311->MapSubwindows();
   fCompositeFrame309->AddFrame(fCanvas311, new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY));

   fRootCanvas284->AddFrame(fCompositeFrame309, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   // status bar
   TGStatusBar *fStatusBar324 = new TGStatusBar(fRootCanvas284,10,16);
   Int_t partsusBar324[] = {33,10,10,47};
   fStatusBar324->SetParts(partsusBar324,4);
   fRootCanvas284->AddFrame(fStatusBar324, new TGLayoutHints(kLHintsLeft | kLHintsBottom | kLHintsExpandX,2,2,1,1));

   fRootCanvas284->SetWindowName("tcnvRane0");
   fRootCanvas284->SetIconName("tcnvRane0");
   fRootCanvas284->SetIconPixmap("macro_s.xpm");
   fRootCanvas284->SetClassHints("ROOT","Canvas");
   fRootCanvas284->SetMWMHints(kMWMDecorAll,
                        kMWMFuncAll,
                        kMWMInputModeless);
   fRootCanvas284->MapSubwindows();
   fHorizontal3DLine302->UnmapWindow();
   fDockableFrame303->UnmapWindow();
   fHorizontal3DLine308->UnmapWindow();
   fCompositeFrame310->UnmapWindow();
   fStatusBar324->UnmapWindow();

   fRootCanvas284->Resize(fRootCanvas284->GetDefaultSize());
   fRootCanvas284->MapWindow();
   fRootCanvas284->Resize(800,536);
}  

import { Routes } from '@angular/router';
import { HomeComponent } from './home/home.component';
import { NewCaseComponent } from './newppkt/newppkt.component';
import { PhenotypeProfileVisualizerComponent } from './visualize/visualize.component';
import { GeneDiseaseComponent } from './util/genedisease/genedisease.component';
import { HelpComponent } from './help/help.component';
import { NotebookComponent } from './notebook/notebook.component';


export const appRoutes: Routes = [
  { path: '', redirectTo: 'home', pathMatch: 'full' },
  { path: 'home', component: HomeComponent },
  { path: 'newppkt', component: NewCaseComponent },
  { path: 'visualize', component: PhenotypeProfileVisualizerComponent },
  { path: 'genedisease', component: GeneDiseaseComponent },
  { path: 'notebook', component: NotebookComponent},
  { path: 'help', component: HelpComponent },
  
  
  // The wildcard safety net stays at the bottom
  { path: '**', redirectTo: 'home' },
];
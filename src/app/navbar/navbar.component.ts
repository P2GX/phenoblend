import { Component, OnInit } from '@angular/core';
import { RouterLink, RouterLinkActive } from '@angular/router';

@Component({
  selector: 'app-navbar',
  standalone: true,
  templateUrl: './navbar.component.html',
  styleUrls: ['./navbar.component.css'],
  imports: [RouterLink, RouterLinkActive],
})
export class NavbarComponent implements OnInit{
  tabs = [
    { label: 'Home', path: '/home' },
    { label: 'Load HPO data', path: '/newppkt' },
    { label: 'Gene/Disease data', path: '/genedisease' },
    { label: 'Visualize', path: '/visualize'},
    { label: 'Notebook', path: '/notebook'}
  ];

  ngOnInit() {
    /* no op */
  }

  isDisabled(tab: { path: string }) {
    console.log("Add is disabled logic if needed");
    return false; 
  }
  logClick(path: string) {
    // no-op console.log(`Link clicked for target path: ${path}`);
  }

}


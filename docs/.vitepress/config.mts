import { defineConfig } from 'vitepress'

export default defineConfig({
  title: "phenoblend",
  description: "Phenoblend Toolkit",
  base: '/phenoblend/', 
  
  themeConfig: {
    nav: [
      { text: 'Home', link: '/' }
    ],
    sidebar: [
      {
        text: 'Getting Started',
        items: [
          { 
            text: 'Setup', 
            items: [
              { text: 'Installation', link: '/setup/installation' },
              { text: 'Developers', link: '/setup/developers' }
            ]
          },
          { 
            text: 'Tutorial', 
            link: '/tutorial/index',
            items: [
              { text: 'HPO(A) Data', link: '/tutorial/data' },
              { text: 'Input data', link: '/tutorial/input' },
              { text: 'Visualizing data', link: '/tutorial/visualizations' },
              { text: 'Overlap plot', link: '/tutorial/overlap' },
              { text: 'Upset plot', link: '/tutorial/upset' },
              { text: 'Phenotypic profile plot', link: '/tutorial/profile' }
            ]
          }
        ]
      }
    ],
    socialLinks: [
      { icon: 'github', link: 'https://github.com/P2GX/phenoblend' }
    ]
  }
})
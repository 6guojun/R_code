import Vue from 'vue'
import Router from 'vue-router'
import Sidebar from '@/components/common/Sidebar'

Vue.use(Router)

export default new Router({
  routes: [
    {
      path: '/',
      name: 'Sidebar',
      component: Sidebar
    }
  ]
})

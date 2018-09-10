<template>
  <div class="container-fluid">
    <div class="row navbar">
      <navbar :navitemLists="navitemLists" :website="website"></navbar>
    </div>
    <div class="row content-container">
        <div class="sidebar" :class="{'col-lg-1 col-md-2 col-sm-3 col-xs-4': !isActive, 'col-lg-2 col-md-3 col-sm-4 col-xs-5': isActive}">
            <sidebar :apiUrl="sidebarApiUrl" :treeviewMode="true"
                     :inputSidebarData="inputSidebarData"
                     :defaultAction=true
                     @activateEvent="setActive"
                     @postLinkEvent="setApiUrl('sidebarApiUrl', $event)">
            </sidebar>
        </div>
        <div v-if="isActive" class="table-grid" 
             :class="{'col-lg-10 col-md-9 col-sm-8 col-xs-7': isActive, 'col-lg-11 col-md-10 col-sm-9 col-xs-8': !isActive}">
            <monitor-plot></monitor-plot>
            <nordata-vue-table :apiUrl="tableApiUrl" :placeholder="placeholder"></nordata-vue-table>
        </div>
    </div>
  </div>
</template>

<script>
import Sidebar from "./components/common/Sidebar"
import Navbar from "./components/common/Navbar"
import NordataVueTable from "./components/hostmgt/NordataVueTable"
import MonitorPlot from "./components/hostmgt/MonitorPlot"

export default {
  name: 'app',
  data () {
    return {
      // sidebarApiUrl: "http://localhost:3000/db",
      placeholder: "UUID",
      sidebarApiUrl: "http://localhost:3000/sidebarData",
      tableApiUrl: "http://localhost:3000/tables",
      isActive: false,  // isActive：激活侧边栏
      website: {
        href: 'http://www.nordata.com.cn',
        name: 'SuperSAN'
      },
      inputSidebarData: null,
      // navbar组件
      navitemLists: {
        project: {
          href: '#project',
          name: '项目管理'
        },
        workflow: {
          href: '#add-workflow',
          name: '新建工作流'
        },
        task: {
          href: '#task',
          name: '查看任务'
        },
        results: {
          href: '#results',
          name: '查看结果'
        },
        runinfo: {
          href: '#runinfo',
          name: '查看监控'
        },
        filemgmt: {
          href: '#filemgmt',
          name: '文件管理'
        },
        help: {
          href: '#help',
          name: '帮助'
        }
      },
    }
  },
  methods: {
    setActive: function(event){
      this.isActive = event
    },
    setApiUrl: function(obj_name, event){
      if(event.link){
        this.sidebarApiUrl = event.link
      }
    }
  },
  components: {
    Sidebar,
    Navbar,
    NordataVueTable,
    MonitorPlot
  }
}
</script>

<style>
html, body {
  height: 100%;
  overflow: hidden;
}

.container-fluid {
  width: 100%;
  height: 100%;
  padding: 0 0;
}

.navbar {
  height: 5%;
  margin-bottom: 0px;
  border: 0px;
}

.content-container {
  height: 95%;
}

.table-grid {
  overflow: scroll;
}

/* just for Chrome */
.table-grid::-webkit-scrollbar {
  display:none
}

.col-lg-1, .col-lg-2, .col-lg-3,
.col-lg-4, .col-lg-5, .col-lg-6,
.col-lg-7, .col-lg-8, .col-lg-9,
.col-lg-10 {
  height: 100%;
  margin-right: -15px;
}
</style>

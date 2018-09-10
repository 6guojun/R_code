<template>
    <aside class="sidebar-container">
        <section  class="sidebar" :class="{'show-menu': showMenu}">
            <ul class="sidebar-menu">
                <!-- Header -->
                <li v-if="sidebarData.header" class="sidebar-header">
                    <a :href="sidebarData.header.link" 
                    @click="openLink(sidebarData.header.link, $event)">{{ sidebarData.header.title }}</a>
                </li>

                <!-- Dashboard -->
                <li v-if="sidebarData.dashboard" class="treeview">
                    <div class="treeview-content">
                        <i class="fa fa-dashboard"></i> 
                        <a :href="sidebarData.dashboard.link"
                        @click="openLink(sidebarData.dashboard.link, $event)"><span>{{ sidebarData.dashboard.title }}</span></a>
                    </div>
                </li>

                <!-- Treeview -->
                <li v-for="treeitem in sidebarData.treeitems" class="treeview">
                    <!-- treeview -->
                    <div class="treeview-content"  v-if="treeviewMode"
                    @click.stop.prevent="openTreeMenu(treeitem.id, $event)">
                    <i class="fa fa-cubes"></i>
                    <a :href="treeitem.link" v-if="treeviewMode" 
                    @click="openLink(treeitem.link, $event)"><span>{{ treeitem.title }}</span></a>
                    <i class="fa fa-angle-left pull-right"
                    :class="{'active': showMenus[treeitem.id]}"></i>
                </div>

                <!-- checkbox -->
                <div class="treeview-content"  v-if="checkboxMode">
                    <input type="checkbox" value="test"/>
                    <a><span>{{ treeitem.title }}</span></a>
                    <i class="fa fa-angle-left pull-right"
                    :class="{'active': showMenus[treeitem.id]}"
                    @click="openTreeMenu(treeitem.id, $event)"></i>
                </div>

                <!-- Secondary Menus -->
                <ul v-if="treeitem.secondmenus" class="sec-treeview-menu"
                :class="{'showtreeview-menu': showMenus[treeitem.id]}">
                <li v-for="treemenuitem in treeitem.secondmenus" class="sec-treeview-menu-item">
                    <!-- treeview -->
                    <div class="sec-treeview-menu-item-content" v-if="treeviewMode"
                    @click.stop.prevent="openTreeMenu(treemenuitem.id, $event)">
                    <a :href="treemenuitem.link" @click="openLink(treemenuitem.link, $event)">
                        <i class="fa fa-circle"></i>
                        <span>{{ treemenuitem.title }}</span>
                    </a>
                    <i class="fa fa-angle-left pull-right"
                    :class="{'active': showMenus[treemenuitem.id]}"></i>
                </div>

                <!-- checkbox -->
                <div class="sec-treeview-menu-item-content" v-if="checkboxMode">
                    <input type="checkbox" value="test" @click="checkBoxes($event)"/>
                    <a><span>{{ treemenuitem.title }}</span></a>
                    <i class="fa fa-angle-left pull-right"
                    :class="{'active': showMenus[treemenuitem.id]}"
                    @click="openTreeMenu(treemenuitem.id, $event)"></i>
                </div>

                <!-- Thirdary Menus -->
                <ul v-if="treemenuitem.thirdmenus" class="third-treeview-menu"
                :class="{'showtreeview-menu': showMenus[treemenuitem.id]}">
                <li v-for="subtreemenuitem in treemenuitem.thirdmenus" class="third-treeview-menu-item">
                    <!-- treeview -->
                    <div class="third-treeview-menu-content" v-if="treeviewMode">
                        <i class="fa fa-circle-o"></i>
                        <a :href="subtreemenuitem.link" v-if="treeviewMode"
                        @click="openLink(subtreemenuitem.link, $event)"><span>{{ subtreemenuitem.title }}</span></a>
                    </div>

                    <!-- checkbox -->
                    <div class="third-treeview-menu-content" v-if="checkboxMode">
                        <input type="checkbox" value="test" @click="checkBoxes($event)"/>
                        <a><span>{{ subtreemenuitem.title }}</span></a>
                    </div>
                </li>
            </ul>
            <!-- Thirdary Menus End -->
        </li>
    </ul>
    <!-- Secondary Menus End -->
</li>
<!-- Treeview -->
</ul>
<button class="close-button" id="close-button" @click="closeMenu"></button>
</section>
<button class="menu-button" id="menu-button" @click="openMenu" :class="{'show-menu-btn': showMenuBtn}"></button>
</aside>
</template>

<script>
    import axios from 'axios'

    export default {
        name: 'sidebar',
        props: {
            apiUrl: {
                type: String,
                required: true
            },
            treeviewMode: Boolean,
            checkboxMode: Boolean,
            eventBus: {
                type: Object,
                required: false
            },
            defaultAction: {
                type: Boolean,
                required: false,
                default: false
            },
            inputSidebarData: null,
        },
        data (){
            return {
                showMenu: false,
                showMenuBtn: true,
                showMenus: [],
                sidebarData: {
                    'header': null,
                    'dashboard': null,
                    'treeitems': null
                }
            }
        },
        methods: {
            openTreeMenu: function(index, event){
                console.log(index)
                for(var i = 1; i <= index; i++){
                    if(this.showMenus[i] === undefined){
                        this.showMenus[i] = false
                    }
                }
                console.log(this.showMenus[index])
                var val = this.showMenus[index]
                this.showMenus.splice(index, 1, !val)
            },
            checkBoxes: function(index, event){

            },
            openMenu: function(){
                this.showMenu = true
                this.showMenuBtn = false
                this.$emit("activateEvent", Boolean(true))
            },
            closeMenu: function(){
                this.showMenu = false
                this.showMenuBtn = true
                this.$emit("activateEvent", Boolean(false))
            },
            openLink: function(link, event){
                event['link'] = link
                if(this.defaultAction){
                    return true
                } else if(this.eventBus){
                    this.eventBus.$emit("postLinkEvent", event)
                    event.preventDefault();
                } else {
                    this.$emit("postLinkEvent", event)
                    event.preventDefault();
                }
            },
            getData: function(){
                console.log("获取Sidebar数据", this.api)
                const vm = this
                axios.get(this.apiUrl)
                .then(res => {
                    console.log(res.data.sidebarData)
                    // vm.sidebarData = res.data.sidebarData
                    vm.sidebarData = res.data
                })
                .catch(function(err) {console.log(err)})
            }
        },
        watch: {
            apiUrl: function(){
                console.log("更新Sidebar数据")
                this.getData()
            }
        },
        created: function(){
            if(this.apiUrl != ""){
                console.log("获取Sidebar数据", this.apiUrl)
                this.getData()
            }else if(this.inputSidebarData){
                this.sidebarData = this.inputSidebarData
            }
        }
    }
</script>

<!-- Add "scoped" attribute to limit CSS to this component only -->
<style scoped>
    .sidebar-container {
        height: 100%;
    }

    .sidebar {
        display: none;
        height: 100%;
        overflow: auto;
        z-index: 810;
        background-color: #222d32;
    }

    .sidebar-menu {
        list-style: none;
        margin: 0;
        padding: 0;
        background-color: #222d32; 
    }

    .sidebar-header, .treeview {
        padding: 12px 5px 12px 15px;
        position: relative;
        margin: 0;
    }

    .sidebar-header {
        font-size: 1.5em;
        color: #4b646f;
        background: #1a2226; 
    }

    .sec-treeview-menu, .sec-treeview-menu-item {
        padding: 12px 0px 0px 5px;
    }

    .third-treeview-menu, .third-treeview-menu-item {
        padding: 12px 0px 0px 5px;
    }

    .sec-treeview-menu, .third-treeview-menu {
        display: none;
    }

    .fa-angle-left {
        width: 20px;
        height: auto;
        padding: 0;
        margin-top: 3px; 
    }

    .active {
        transform: rotate(-90deg); 
        padding-bottom: 10px;
        padding-left: 10px;
    }

    .treeview-content, 
    .sec-treeview-menu-item-content,
    .third-treeview-menu-content {
        cursor: pointer;
        text-decoration: none;
        border-left: 3px solid transparent;
        color: #b8c7ce;
        padding: 0 0;
        display: block;
        border-left: 0;
        font-size: 1.2em;
    }

    a {
        cursor: pointer;
        text-decoration: none;
        color: #b8c7ce;
    }

    a:hover {
        color: #fff;
        border-left-color: #3c8dbc; 
    }

    /* Menu Button */
    .menu-button {
        position: fixed;
        z-index: 1000;
        margin: 0.3em 1em;
        padding: 0;
        width: 2.2em;
        height: 2em;
        border: none;
        text-indent: 2.5em;
        font-size: 1.5em;
        color: transparent;
        background: transparent;
        outline: none;
        display: none;
    }

    .menu-button::before {
        position: absolute;
        top: 0.5em;
        right: 0.5em;
        bottom: 0.5em;
        left: 0em;
        background: linear-gradient(#373a47 20%, transparent 20%, 
            transparent 40%, #373a47 40%, 
            #373a47 60%, transparent 60%, 
            transparent 80%, #373a47 80%);
        content: '';
    }

    .menu-button:hover {
        opacity: 0.6;
    }

    /* Close Button */
    .close-button {
        width: 1em;
        height: 1em;
        position: absolute;
        right: 2em;
        top: 1em;
        overflow: hidden;
        text-indent: 1em;
        font-size: 0.75em;
        border: none;
        background: transparent;
        color: transparent;
    }

    .close-button::before,
    .close-button::after {
        content: '';
        position: absolute;
        width: 3px;
        height: 100%;
        top: 0;
        left: 50%;
        background: #bdc3c7;
    }

    .close-button::before {
        -webkit-transform: rotate(45deg);
        transform: rotate(45deg);
    }

    .close-button::after {
        -webkit-transform: rotate(-45deg);
        transform: rotate(-45deg);
    }

    /* Shown menu */
    .show-menu .sidebar {
        -webkit-transform: translate3d(0,0,0);
        transform: translate3d(0,0,0);
        -webkit-transition: -webkit-transform 0.8s;
        transition: transform 0.8s;
        -webkit-transition-timing-function: cubic-bezier(0.7,0,0.3,1);
        transition-timing-function: cubic-bezier(0.7,0,0.3,1);
    }

    .show-menu {
        display: block;
    }

    .show-menu-btn {
        display: block;
    }

    .showtreeview-menu {
        display: block;
        list-style: none;
    }
</style>

<template>
    <tabs :tabsData=tabsData>
        <button type='button' class='close' data-dismiss='alert' 
        @click="closeTabs" slot="button">x</button>
        <tab v-for="(tabContent, tabsKey) in tabsData" :key="tabContent.id"
        :name="tabContent.name" :title="tabContent.title"
        :active="tabContent.active.toString() == 'true'">
        <p v-html='tabContent.content'></p>
    </tab>
</tabs>
</template>

<script>
// CSS Style
import 'bootstrap/dist/css/bootstrap.css'
import 'bootstrap/dist/css/bootstrap-theme.css'

// External Components
import axios from 'axios'
import Vue from "vue"

// Custom Components - Common
import Tabs from '../common/Tabs'
import Tab from '../common/Tab'

export default {
    components: {
        Tabs,
        Tab
    },
    props: {
        apiUrl: {
            type: String,
            required: false
        },
        inputTabsData: null
    },
    data(){
        return {
            tabsData: null
        }
    },
    methods: {
        getData: function(){
            const vm = this
            axios.get(this.apiUrl)
            .then(res => {
                vm.tabsData = res.data.tabsData
                console.log(res.data.tabsData)
            })
            .catch(function(err) {console.log(err)})
        },
        closeTabs(event) {
            this.$events.fire('close-tabs')
        }
    },
    created: function(){
        if(this.apiUrl != ""){
            this.getData()
        }else if(this.inputTabsData){
            this.tabsData = this.inputTabsData
        }
    }
}
</script>

<style scoped>
    .tabs-close {
        display: none;
    }
</style>
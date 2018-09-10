<template>
    <div class="container">
        <slot name="button"></slot>
        <ul class="nav nav-tabs">
            <li role="presentation" :key="tabContent.id"
            :class="{'active': tabContent.active.toString() == 'true'}"
            v-for="(tabContent, tabKey) in tabsData">
            <a @click="showTabContent"
            :aria-controls="tabContent.name" 
            role="tab" data-toggle="tab" :name="tabKey">
            {{tabContent.title}}
        </a>
    </li>
</ul>
<div class="tab-content">
    <slot> </slot>
</div>
</div>
</template>
<script>
    export default{
        props: {
            tabsData: {
                type: Object
            }
        },
        methods: {
            showTabContent(event) {
                console.log('showTabContent:', event.target.name)
                console.log(this.tabsData[event.target.name])
                for(var key in this.tabsData){
                    console.log(key)
                    if(key != event.target.name){
                        this.tabsData[key]['active'] = false
                    }else{
                        this.tabsData[key]['active'] = true
                    }
                }
            }
        }
    }
</script>
<style>
    .nav-tabs {
        background-color: #eeeeee;
        padding-top:5px;
        padding-left:5px;
        border: 1px solid #ddd;

    }

    .tab-content {
        color : #000000;
        background-color: #ffffff;
        padding : 5px 15px;
        margin-top: -1px;
        border: 1px solid #ddd;
    }

    button.close {
        padding-right: 5px;
    }

    .nav-tabs li:hover {
        cursor: pointer;
    }
</style>

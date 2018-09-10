<template>
    <div class="ui container">
        <filter-bar :placeholder="placeholder" :label="filterBarLabel"
                    :searchBtnLabel="searchBtnLabel"
                    :resetBtnLabel="resetBtnLabel"
                    class="float-right">
        </filter-bar>
        <vuetable-pagination-dropdown class="float-left">
        </vuetable-pagination-dropdown>
        <vuetable ref="vuetable" :api-url="apiUrl" class="vuetable"
                  :fields="fields" :append-params="moreParams"
                  :per-page="1" @vuetable:cell-clicked="onCellClicked"
                  pagination-path="" 
                  @vuetable:pagination-data="onPaginationData">
        </vuetable>
        <div class="vuetable-pagination ui basic segment grid">
            <vuetable-pagination-info ref="paginationInfoBottom" 
                                      :infoTemplate="infoTemplate">
            </vuetable-pagination-info>
            <vuetable-pagination ref="paginationBottom"
                                 @vuetable-pagination:change-page="onChangePage">
            </vuetable-pagination>
        </div>
        <div class="tabs-container" v-if="isActive">
            <nordata-vue-tabs class="content-container" :api-url="tabsApiUrl">
            </nordata-vue-tabs>
            <div class="mask-window" @click="onCloseTabs"></div>
        </div>
    </div>
</template>

<script type="text/javascript">
    // CSS Stype
    import '../../assets/css/semantic-ui.css'
    import 'bootstrap/dist/css/bootstrap.css'
    import 'bootstrap/dist/css/bootstrap-theme.css'

    // External Components
    import Vue from "vue"
    import VueEvents from 'vue-events'

    // Custom Components - Common
    import FilterBar from '../common/FilterBar'
    import Tab from '../common/Tab'
    import Vuetable from "../common/Vuetable"
    import VuetablePagination from "../common/VuetablePagination"
    import VuetablePaginationInfo from "../common/VuetablePaginationInfo"
    import VuetablePaginationDropdown from "../common/VuetablePaginationDropdown"

    // Custom Components
    import NordataVueTabs from './NordataVueTabs'
    import CustomActions from "./CustomActions"

    Vue.use(VueEvents)
    Vue.component("custom-actions", CustomActions)

    export default {
        components: {
            Vuetable,
            FilterBar,
            NordataVueTabs,
            VuetablePagination,
            VuetablePaginationInfo,
            VuetablePaginationDropdown,
        },
        props: {
            apiUrl: {
                type: String,
                required: true
            },
            placeholder: {
                default: "UUID/主机名/组名",
                type: String,
                required: false
            }
        },
        data() {
            return {
                fields: [
                    {
                        name: "__sequence",
                        title: "#",
                        titleClass: "center aligned",
                        dataClass: "right aligned"
                    },
                    {
                        name: "powerstatus",
                        title: "电源状态",
                        sortField: "powerstatus",
                        callback: "labelStatus",
                        titleClass: "center aligned",
                        dataClass: "center aligned"
                    },
                    {
                        name: "hostuuid",
                        title: "UUID",
                        sortField: 'hostuuid'
                    },
                    {
                        name: "hostname", 
                        title: "主机名",
                        sortField: "hostname"
                    },
                    {
                        name: "ipmi_ipaddr",
                        title: "IPMI IP地址",
                        sortField: "ipmi_ipaddr"
                    },
                    {
                        name: "ipmi_mac",
                        title: "IPMI MAC地址",
                        sortField: 'ipmi_mac'
                    },
                    {
                        name: "group_name",
                        title: "组名",
                        callback: "groupLabel",
                        sortField: "group_name"
                    },
                    {
                        name: "__component:custom-actions",
                        title: "动作",
                        titleClass: "center aligned",
                        dataClass: "center aligned"
                    }
                ],
                moreParams: {},
                tabsApiUrl: "http://localhost:3000/tabs",
                isActive: false,
                infoTemplate: "页数：{from}/{to} 页 总数：{total}条",
                filterBarLabel: "过滤器",
                searchBtnLabel: "搜索",
                resetBtnLabel: "重置"
            }
        },
        mounted() {
            this.$events.$on('filter-set', eventData => this.onFilterSet(eventData))
            this.$events.$on('filter-reset', e => this.onFilterReSet())
            this.$events.$on('close-tabs', e => this.onCloseTabs())
            this.$events.$on('show-tabs', eventData => this.onShowTabs(eventData))
        },
        methods: {
            onCloseTabs(rowData){
                console.log(this.isActive, "rowData", rowData)
                this.isActive = false
            },
            onShowTabs(event){
                console.log(this.isActive)
                this.isActive = true
            },
            onCellClicked(rowData, cell, event) {
                console.log('cellClicked: ', cell, rowData)
                console.log(this.isActive)
                this.isActive = true
            },
            onFilterSet(filterText) {
                console.log('filter-set', filterText)
                this.moreParams = {
                    'filter': filterText
                }
                Vue.nextTick( () => this.$refs.vuetable.refresh())
            },
            onFilterReSet() {
                console.log('filter-reset')
                this.moreParams = {}
                Vue.nextTick( () => this.$refs.vuetable.refresh())
            },
            groupLabel(value) {
                return value == "management"
                ? "<span class='label label-info'><i class='glyphicon glyphicon-star'></i> 管理组 </span>"
                : "<span class='label label-danger'><i class='glyphicon glyphicon-heart'></i> 计算组 </span>"
            },
            labelStatus(value) {
                var newValue = value.toUpperCase()
                if (newValue == "POWER_ON") {
                    return "<span class='label label-success'>开机</span>"
                } else if (newValue == "POWER_OFF") {
                    return "<span class='label label-primary'>关机</span>"
                } else if (newValue == "REBOOT") {
                    return "<span class='label label-warning'>重启</span>"
                } else {
                    return "<span class='label label-danger'>未知</span>"
                }
            },
            onPaginationData(paginationData) {
                this.$refs.paginationBottom.setPaginationData(paginationData)
                this.$refs.paginationInfoBottom.setPaginationData(paginationData)
            },
            onChangePage(page) {
                this.$refs.vuetable.changePage(page)
            }
        }
    }
</script>

<style scoped>
    .tabs-container {
        display: flex;
        position: fixed;
        top: 20%;
        left: 20%;
        align-items:center;
        width: 60%;
        height: 60%;
        justify-content: center;
    }
    .mask-window {
        z-index: 999;
        position: fixed;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background: #000;
        opacity: 0.2;
    } 
    .vuetable {
        border-radius: 5px!important;
    }
    .ui.container {
        width: 100%!important;
        margin-bottom: 15px;
    }
    .content-container {
        z-index: 1000;
    }
    .float-right {
        position: relative;
        float: right!important;
    }
    .float-left {
        position: relative;
        float: left!important;
    }
    .large.text {
       font-size: 1.2rem;
    }
</style>
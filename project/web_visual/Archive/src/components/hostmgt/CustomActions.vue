<template>
    <div class="custom-actions">
        <button class="btn btn-default" @click="itemAction('power-on-item', rowData, rowIndex)">开机</button>
        <button class="btn btn-default" @click="itemAction('power-off-item', rowData, rowIndex)">关机</i></button>
        <button class="btn btn-default" @click="itemAction('reboot-item', rowData, rowIndex)">重启</button>
        <button class="btn btn-default" @click="itemAction('details-item', rowData, rowIndex)">详细信息</button>
        <button class="btn btn-default" @click="itemAction('monitor-item', rowData, rowIndex)">监控信息</button>
    </div>
</template>

<script type="text/javascript">
    // CSS Style
    import 'bootstrap/dist/css/bootstrap.css'
    import 'bootstrap/dist/css/bootstrap-theme.css'

    export default {
        components: {

        },
        name: 'custom-actions',
        props: {
            rowData: {
                type: Object,
                required: true
            },
            rowIndex: {
                type: Number
            }
        },
        methods: {
            itemAction(action, rowData, index) {
                console.log(action, rowData, index)
                if(action == "details-item"){
                    this.$events.fire("show-tabs", rowData)
                }
            }
        }
    }
</script>

<style scoped>
    .custom-actions button.ui.button {
        padding: 8px 8px;
    }
    .custom-actions button.ui.button > i.icon {
        margin: auto !important;
    }
</style>